import dxpy
import inspect
import os
import re
from threading import Lock
from typing import Optional
from dxpy.api import project_describe

from dxcint.Logger import Logger


class ContextError(Exception):
    """
    Class to handle Context errors
    """


class Context(object):
    def __init__(
        self,
        project: str,
        repo_root: str,
        folder: Optional[str] = None,
        logger_verbosity: str = "info",
    ):
        """
        Immutable class maintaining the passed arguments and values which dxcint was initiated with
        Args:
            project: project on the platform where the tests will be run
            repo_root: directory path of the dxCompiler root
            folder: directory name on the platform to build dxCompiler assets. If None: /builds/<USER_NAME>/<DXC_VERISON>
            logger_verbosity: verbosity level for Logger class. See dxcint.Logger
        """
        self._logger = Logger.make(name=__name__, verbosity=logger_verbosity)
        self._project_id = self._resolve_project(project=project)
        self._user = dxpy.whoami()
        self._repo_root_dir = os.path.realpath(repo_root)
        self._compiler_version = self._get_version()
        self._platform_build_dir = folder = (
            folder or f"/builds/{self._user}/{self._compiler_version}"
        )
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

    @property
    def logger(self):
        return self._logger

    @staticmethod
    def _resolve_project(project: str) -> str:
        # First, see if the project is a project-id.
        try:
            project_obj = dxpy.DXProject(project)
            return project_obj.get_id()
        except dxpy.DXError:
            found_projects = list(dxpy.find_projects(name=project, level="CONTRIBUTE"))
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


class ContextEmpty(Context):
    def __init__(self):
        self._project_id = ""
        self._user = None
        self._repo_root_dir = ""
        self._platform_build_dir = "."
        self._compiler_version = ""
        self._lock = Lock()
        self._project_info = None
        self._logger = Logger.make(name=__name__, verbosity="error")

    def __setattr__(self, *args):
        if inspect.stack()[1][3] == "__init__":
            object.__setattr__(self, *args)
        else:
            raise ContextError("Context class is immutable")
