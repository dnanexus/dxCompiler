import os
import re
import json
import tarfile
import subprocess as sp

from typing import Dict
from pathlib import Path

from dxcint.Context import Context


class DependencyError(Exception):
    """
    Class to handle Dependency exceptions
    """


class DependencyFactory(object):
    def __init__(self, config_file: str, context: Context):
        self._type_switch = {
            "binary": BinaryDependency,
            "tarball": TarballDependency,
            "package_manager": PackageDependency,
            "default": Dependency,
        }
        self._config_file = config_file
        self._context = context

    def make(self):
        detected_type = self._detect_type()
        return self._type_switch[detected_type](self._config_file, self._context)

    def _detect_type(self) -> str:
        """
        This logic defines the precedence of detection. First precedence: If config has `package_manager` field, then
        it's a package dependency. Later this will be reconsidered to always deliver dependencies as packages

        Returns: str. Type of package.
        """
        with open(self._config_file, "r") as config_handle:
            config = json.load(config_handle)
        if config.get("package_manager", None):
            return "package_manager"

        src_link = config.get("source_link", None)
        if src_link:
            if re.search(r"tar\..z2?$", src_link):
                return "tarball"
            else:
                return "binary"
        else:
            return "default"


class Dependency(object):
    def __init__(self, config_file: str, context: Context):
        """
        Parent class for providing the necessary dependencies for building dxCompiler assets on the platform
        Args:
            config_file: config file for a dependency. Should be located in ${DXCOMPILER_REPO}/dxcint/dependencies/config
            context: instance of dxcint.Context
        """
        self._dependencies_root = Path(os.path.join(config_file, "../..")).resolve()
        with open(config_file, "r") as config_handle:
            config = json.load(config_handle)
        # Required config parameters
        self._name = config.get("name")
        self._languages = config.get("languages")
        self._version = config.get("version")
        self._source_link = config.get("source_link", None)
        self._package_manager = config.get("package_manager", None)
        self._env = config.get("env", None)
        self._context = context

    @property
    def languages(self):
        return self._languages

    def link(self, symlink_destination: Path) -> None:
        pass

    def export_spec(self) -> None:
        pass


class BinaryDependency(Dependency):
    def __init__(self, config_file: str, context: Context):
        super().__init__(config_file, context)
        # additional attributes
        self._local_dir = self._get_exec()
        if not self._source_link:
            raise DependencyError(
                "BinaryDependency(): `source_link` field in the config file can not be None for BinaryDependency class"
            )

    def link(self, link_destination_dir: Path) -> Path:
        """
        Method to create a symlink for the dependency executable in the asset resources
        Args:
            link_destination_dir: Path. Path to the link.

        Returns: Path. Path of a symlink

        """
        if not os.path.exists(link_destination_dir):
            os.makedirs(link_destination_dir)
        link_path = Path(os.path.join(link_destination_dir, Path(self._local_dir).name))
        if not os.path.exists(link_path):
            os.link(self._local_dir, link_path)
        return link_path

    def _get_exec(self) -> Path:
        """
        Method to download (if not available) the dependency executable

        Returns: Path. Location of a downloaded binary exec
        """
        local_dependency = os.path.join(
            self._dependencies_root, self._name, self._version, self._name
        )
        if os.path.exists(local_dependency):
            self._context.logger.info(
                f"BinaryDependency._get_exec(): Found binary in local directory {local_dependency}. "
                f"Will be using that one."
            )
            return Path(local_dependency)
        else:
            os.makedirs(os.path.dirname(local_dependency), exist_ok=True)
            self._context.logger.info(
                f"BinaryDependency._get_exec(): Downloading {self._name} version {self._version}"
            )
            download_cmd = ["wget", self._source_link, "-O", local_dependency]
            proc = sp.Popen(download_cmd, stdout=sp.PIPE, stderr=sp.PIPE)
            out, err = proc.communicate()
            if proc.returncode != 0:
                raise DependencyError(f"Executable download raised {err.decode()}")
            else:
                self._context.logger.info(
                    f"BinaryDependency._get_exec(): Executable download returned {out.decode()}\n{err.decode()}"
                )
            os.chmod(local_dependency, 0o775)
        return Path(local_dependency)


class TarballDependency(BinaryDependency):
    def __init__(self, config_file: str, context: Context):
        super().__init__(config_file, context)
        # additional attributes
        self._local_dir = self._get_exec()
        if not self._source_link:
            raise DependencyError(
                "TarballDependency(): `source_link` field in the config file can not be None for "
                "TarballDependency class"
            )

    def _get_exec(self) -> Path:
        """
        Method to download (if not available) and unpack the dependency executable

        Returns: Path. Location of an extracted binary exec
        """
        local_dependency_tarball = os.path.join(
            self._dependencies_root, self._name, self._version, self._name
        )
        if not os.path.exists(local_dependency_tarball):
            os.makedirs(os.path.dirname(local_dependency_tarball), exist_ok=True)
            self._context.logger.info(
                f"Downloading tarball for {self._name} version {self._version}"
            )
            download_cmd = ["wget", self._source_link, "-O", local_dependency_tarball]
            proc = sp.Popen(download_cmd, stdout=sp.PIPE, stderr=sp.PIPE)
            out, err = proc.communicate()
            if proc.returncode != 0:
                raise DependencyError(f"Executable download raised {err.decode()}")
            else:
                self._context.logger.info(
                    f"Executable download returned {out.decode()}\n{err.decode()}"
                )
        else:
            self._context.logger.info(
                f"TarballDependency._get_exec(): Found tarball in local directory {local_dependency_tarball}. "
                f"Will be using that one."
            )
        with tarfile.open(local_dependency_tarball, "r") as tar_handle:
            top_tar_dir = tar_handle.getnames()[0]
            if not os.path.exists(
                os.path.join(os.path.dirname(local_dependency_tarball), top_tar_dir)
            ):
                tar_handle.extractall(os.path.dirname(local_dependency_tarball))
        local_dependency = os.path.join(
            os.path.dirname(local_dependency_tarball), top_tar_dir, "bin", self._name
        )
        os.chmod(local_dependency, 0o775)
        return Path(local_dependency)


class PackageDependency(Dependency):
    def __init__(self, config_file: str, context: Context):
        super().__init__(config_file, context)
        if not self._package_manager:
            raise DependencyError(
                "PackageDependency(): `package_manager` field in the config file can not be None for "
                "PackageDependency class"
            )

    def export_spec(self) -> Dict:
        return {
            "name": self._name,
            "package_manager": self._package_manager,
            "version": self._version,
        }
