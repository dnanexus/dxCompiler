import logging
import dxpy
import inspect
import os
import re
from threading import Lock
from typing import Optional
from dxpy.api import project_new_folder, project_describe, project_remove_folder


class ContextError(Exception):
    """
    Class to handle Context errors
    """


class Context(object):
    def __init__(
            self,
            project: str,
            repo_root: str,
            folder: Optional[str] = None
    ):
        self._project_id = self._resolve_project(project=project)
        self._user = dxpy.whoami()
        self._repo_root_dir = repo_root
        self._compiler_version = self._get_version()
        self._platform_build_dir = self._create_platform_build_folder(folder)
        self._lock = Lock()
        self._project_info = project_describe(self._project_id)

    def __setattr__(self, *args):
        if inspect.stack()[1][3] == "__init__":
            object.__setattr__(self, *args)
        else:
            raise ContextError("Context class is immutable")

    @property
    def project_id(self):
        return self._project_id

    @property
    def repo_root_dir(self):
        return self._repo_root_dir

    @property
    def platform_build_dir(self):
        return self._platform_build_dir

    @property
    def version(self):
        return self._compiler_version

    @property
    def lock(self):
        return self._lock

    @property
    def project_info(self):
        return self._project_info

    @staticmethod
    def _resolve_project(project: str) -> str:
        # First, see if the project is a project-id.
        try:
            project = dxpy.DXProject(project)
            return project.get_id()
        except dxpy.DXError:
            found_projects = list(
                dxpy.find_projects(name=project, level="VIEW")
            )
            if len(found_projects) == 0:
                raise ContextError(
                    f"Context._resolve_project(): Could not find project `{project}`."
                )
            elif len(found_projects) > 1:
                raise ContextError(
                    f"Context._resolve_project(): found multiple projects with name `{project}`"
                )
            else:
                return found_projects[0].get("id")

    def _get_version(self) -> str:
        application_conf_path = os.path.join(
            self._repo_root_dir, "core", "src", "main", "resources", "application.conf"
        )
        pattern = re.compile(r"^(\s*)(version)(\s*)(=)(\s*)(\S+)(\s*)$")
        with open(application_conf_path, "r") as fd:
            for line in fd:
                line_clean = line.replace('"', "").replace("'", "")
                m = re.match(pattern, line_clean)
                if m is not None:
                    return m.group(6).strip()
        raise Exception(f"version ID not found in {application_conf_path}")

    def _create_platform_build_folder(self, folder: str) -> str:
        folder = folder or f"/builds/{self._user}/{self._compiler_version}"
        _ = self._clean_up_build_dir(folder)
        _ = self._create_build_subdirs(base_dir=folder)
        logging.info(f"Project: {self._project_id}.\nBuild directory {folder}")
        return folder

    def _create_build_subdirs(self, base_dir: str) -> bool:
        """
        Impure function. Creates destination subdirs on the platform
        Args:
            base_dir: str. Parent dir
        Returns: bool. True when executing without errors
        """
        subdirectories = (os.path.join(base_dir, x) for x in ("applets", "test"))
        for subdir in subdirectories:
            _ = project_new_folder(
                object_id=self._project_id,
                input_params={"folder": subdir, "parents": True}
            )
        return True

    def _clean_up_build_dir(self, folder) -> bool:
        """
        Impure function. Cleans up (removes) target build directories on the platform
        Args:
            folder: str. Dir on the platform to clean up.
        Returns: bool. True when executing without errors
        """
        _ = project_remove_folder(
            object_id=self._project_id,
            input_params={"folder": folder, "force": True, "recurse": True}
        )
        return True


class ContextEmpty(Context):
    def __init__(self):
        self._project_id = None
        self._user = None
        self._repo_root_dir = "."
        self._platform_build_dir = None
        self._compiler_version = None
        self._lock = Lock()
        self._project_info = None

    def __setattr__(self, *args):
        if inspect.stack()[1][3] == "__init__":
            object.__setattr__(self, *args)
        else:
            raise ContextError("Context class is immutable")