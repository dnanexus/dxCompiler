import json
import os
import logging
import shutil
from glob import glob
from pathlib import Path
from typing import List, Dict

from dxcint.RegisteredTest import RegisteredTest


class TestDiscoveryError(Exception):
    """
    Class to handle TestDiscovery errors
    """


class TestDiscovery(object):
    def __init__(self, **test_kwargs):
        """
        Class to handle discovery and addition of the tests to the suite
        :param test_kwargs: use only for unit testing of this class
        """
        self._config_location = test_kwargs.get("config", self._resolve_from_root("config"))
        self._resources_location = test_kwargs.get("resources", self._resolve_from_root("resources"))
        self._allowed_wf_extensions = {
            "cwl": 1,
            "cwl.json": 2,
            "wdl": 3
        }

    def discover(self) -> List[RegisteredTest]:
        pass

    def add_tests(self, dir_name: str, extension: str, suite: str, category: str) -> List[str]:
        """
        Adds tests to the test suite.
        :param dir_name: str. Location of the test-related files (workflow, inputs-outputs, etc.)
        :param extension: str. Test workflow format. Allowed are 'cwl', 'cwl.json', 'wdl'
        :param suite: str. Test suite name. Usually a team-defined group of tests to be run in each CI/CD step
        :param category: str. Test category name. Usually a team-defined category that reflects a type of the test
        :return: List[str]. Names of the added tests
        """
        added_tests = []
        if extension not in self._allowed_wf_extensions.keys():
            raise TestDiscoveryError(
                f"TestDiscovery.add_tests(): {extension} is not recognized. "
                f"Allowed extensions are {self._allowed_wf_extensions.keys()}"
            )
        for test_file in glob(os.path.join(dir_name, f"*.{extension}")):
            test_name_parts = os.path.basename(test_file).split(f".{extension}")
            if test_name_parts[1]:
                raise TestDiscoveryError(
                    f"TestDiscovery.add_tests(): test file name {test_file} contains extension in the middle of the"
                    f"string. Please rename your test file name such that extension `{extension}` is only in the end "
                    f"of the file name"
                )
            test_name = test_name_parts[0]
            test_file_with_dependencies = glob(os.path.join(dir_name, f"{test_name}*"))
            self._copy_resource_files(test_file_with_dependencies, category)
            # Next line inefficient I/O, but an attempt to make EVERY test addition atomic
            self._add_name_to_config(test_name, suite, category)
            added_tests.append(test_name)
        logging.warning(f"Added {len(added_tests)} tests to the suite `{suite}` under category `{category}`")
        return added_tests

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

    @staticmethod
    def _resolve_from_root(dir_name: str) -> Path:
        return Path(os.path.join(os.path.dirname(os.path.abspath(__file__)), f"../../{dir_name}")).resolve()

    def _add_name_to_config(self, test_name: str, suite: str, category: str) -> None:
        suite_file_path = os.path.join(self._config_location, f"{suite}.json")
        if not os.path.exists(suite_file_path):
            logging.warning(f"Test suite `{suite}` does not exist. Creating...")
        with open(suite_file_path, "r+") as suite_config_handle:
            suite_config = json.load(suite_config_handle)
            existing_category = suite_config.get(category, test_name)
            extended_test_collection = set.union(set(existing_category), {test_name})
            suite_config.update({category: list(extended_test_collection)})
            suite_config_handle.seek(0)
            json.dump(suite_config, suite_config_handle)
            suite_config_handle.truncate()

    def _copy_resource_files(self, files: List[str], category: str) -> None:
        category_location = os.path.join(self._resources_location, category)
        if not os.path.exists(category_location):
            logging.warning(f"Category location {category_location} does not exist. Creating...")
            os.mkdir(category_location)
        for one_test_file in files:
            destination_file = os.path.join(category_location, os.path.basename(one_test_file))
            if os.path.exists(destination_file):
                logging.warning(f"Test file {one_test_file} exists. Replacing with the new file")
            shutil.copyfile(one_test_file, destination_file)
