import os
from pathlib import Path


class Dependency(object):
    def __init__(self, destination: Path):
        self._dependency_exec = self._get_exec(destination, version)

    def link(self, symlink_destination: Path) -> str:
        """
        Method to create a symlink for the dependency executable in the asset resources
        :param symlink_destination: Path. Dir name of the
        :return: Path. Path of a symlink
        """
        os.link(res, os.path.join(self._local_asset_dirs.get("bin"), Path(res).name))

    def _get_exec(self, version: str):
        pass
