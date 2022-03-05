import json
from typing import List, Dict
from dxcint.RegisteredTest import RegisteredTest


class TestDiscoveryError(Exception):
    """
    Class to handle TestDiscovery errors
    """


class TestDiscovery(object):
    def __init__(self, location: str):
        self._location = location
        self._config = self._read_config()

    def discover(self) -> List[RegisteredTest]:

        pass

    def _read_config(self) -> Dict:
        with open(self._location, "r") as config_handle:
            config = json.load(config_handle)
            if self._config_linter(config):
                return config

    @staticmethod
    def _config_linter(config) -> bool:
        for test_category, tests in config.items():
            if not isinstance(tests, List):
                raise TestDiscoveryError(
                    f"Unexpected config structure for category {test_category}. "
                    f"Expected an array of test names, received {tests}"
                )
            else:
                continue
        return True


