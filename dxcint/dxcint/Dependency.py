import os
from pathlib import Path


class Dependency(object):
    def __init__(self, ):
        self._dependency_resources = self._get_exec(source, version)

    def link(self, symlink_destination: Path, language: str) -> str:
        """
        Method to create a symlink for the dependency executable in the asset resources
        :param symlink_destination: Path. Dir name of the
        :return: Path. Path of a symlink
        """
        link_path = os.path.join(symlink_destination, Path(self._dependency_resources).name)
        os.link(self._dependency_resources, link_path)
        return link_path

    def _get_exec(self, version: str):
        """
        Method to download (if not available) the dependency executable
        :param version:
        :return:
        """
        pass

    def _update_dot_env(self):
        """
        Method to update .env file with current dependency
        :return:
        """
        pass
        dot_env = "\n".join(f"{key}={val}" for key, val in env_vars.items())
        dot_env_file = os.path.join(self._local_asset_dirs.get("home"), ".env")
        with open(dot_env_file, "wt") as dot_env_handle:
            dot_env_handle.write(dot_env)

