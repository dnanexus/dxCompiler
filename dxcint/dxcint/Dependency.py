import logging
import os
import json
import subprocess as sp

from typing import Dict
from pathlib import Path


class DependencyError(Exception):
    """
    Class to handle Dependency exceptions
    """


class DependencyFactory(object):
    def __init__(self, config_file: str):
        self._type_switch = {
            "source": BinaryDependency,
            "package_manager": PackageDependency
        }
        self._config_file = config_file

    def make(self):
        detected_type = self._detect_type()
        return self._type_switch[detected_type](self._config_file)

    def _detect_type(self) -> str:
        """
        This logic defines the precedence of detection. First precedence: If config has `package_manager` field, then
        it's a package dependency. Later this will be reconsidered to always deliver dependencies as packages
        :return:
        """
        with open(self._config_file, "r") as config_handle:
            config = json.load(config_handle)
        if config.get("package_manager"):
            return "package_manager"
        if config.get("source_link"):
            return "source"


class Dependency(object):
    def __init__(self, config_file: str):
        self._dependencies_root = Path(os.path.join(config_file, "../..")).resolve()
        with open(config_file, "r") as config_handle:
            config = json.load(config_handle)
        # Required config parameters
        self._name = config.get("name")
        self._languages = config.get("languages")
        self._version = config.get("version")
        self._source_link = config.get("source_link", None)
        self._package_manager = config.get("package_manager", None)

    @property
    def languages(self):
        return self._languages

    def update_dot_env(self, home_dir: Path) -> bool:
        """
        Method to update .env file with current dependency
        :param home_dir: Path. Perspective "home" directory in asset "resources"
        :return:
        """
        if not str(home_dir).endswith("home/dnanexus"):
            raise DependencyError(
                f"Dependency._update_dot_env(): provided home dir does not end with `home/dnanexus` - .env file will "
                f"not have a desired effect on the platform. Provided path `{str(home_dir)}`"
            )
        dot_env = f"{self._name}={self._version}\n"
        dot_env_file = os.path.join(home_dir, ".env")
        with open(dot_env_file, "a") as dot_env_handle:
            dot_env_handle.write(dot_env)
        return True

    def link(self, symlink_destination: Path) -> None:
        pass

    def export_spec(self) -> None:
        pass


class BinaryDependency(Dependency):
    def __init__(self, config_file: str):
        super().__init__(config_file)
        # additional attributes
        self._local_dir = self._get_exec()
        if not self._source_link:
            raise DependencyError(
                f"BinaryDependency(): `source_link` field in the config file can not be None for BinaryDependency class"
            )

    def link(self, symlink_destination_dir: Path) -> Path:
        """
        Method to create a symlink for the dependency executable in the asset resources
        :param symlink_destination_dir: Path. Dir name of the
        :return: Path. Path of a symlink
        """
        if not os.path.exists(symlink_destination_dir):
            os.makedirs(symlink_destination_dir)
        link_path = Path(os.path.join(symlink_destination_dir, Path(self._local_dir).name))
        os.link(self._local_dir, link_path)
        return link_path

    def _get_exec(self) -> Path:
        """
        Method to download (if not available) the dependency executable
        :return: Path. Location of a downloaded binary exec
        """
        local_dependency = os.path.join(
            self._dependencies_root, self._name, self._version, self._name
        )
        if os.path.exists(local_dependency):
            logging.info(f"Found binary in local directory {local_dependency}. Will be using that one.")
            return Path(local_dependency)
        else:
            os.makedirs(os.path.dirname(local_dependency), exist_ok=True)
            logging.info(f"Downloading {self._name} version {self._version}")
            download_cmd = ["wget", self._source_link, "-O", local_dependency]
            proc = sp.Popen(download_cmd, stdout=sp.PIPE, stderr=sp.PIPE)
            out, err = proc.communicate()
            if proc.returncode != 0:
                raise DependencyError(f"Executable download raised {err.decode()}")
            else:
                logging.info(f"Executable download returned {out.decode()}\n{err.decode()}")
            os.chmod(local_dependency, 0o775)
        return Path(local_dependency)


class PackageDependency(Dependency):
    def __init__(self, config_file: str):
        super().__init__(config_file)
        if not self._package_manager:
            raise DependencyError(
                f"PackageDependency(): `package_manager` field in the config file can not be None for "
                f"PackageDependency class"
            )

    def export_spec(self) -> Dict:
        return {
            "name": self._name,
            "package_manager": self._package_manager,
            "version": self._version
        }
