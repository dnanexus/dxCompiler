import json
import os
import logging
import shutil
from glob import glob
from pathlib import Path
from typing import List, Dict, Tuple, Set

from dxcint.RegisteredTest import RegisteredTest, RegisteredTestFactory
from dxcint.Dependency import DependencyFactory, Dependency
from dxcint.Context import Context
from dxcint.utils import rm_suffix


class TestDiscoveryError(Exception):
    """
    Class to handle TestDiscovery errors
    """


class TestDiscovery(object):
    def __init__(self, context: Context, **test_kwargs):
        """
        Class to handle discovery and addition of the tests to the suite
        Args:
            context: Context. Global context of the suite
            **test_kwargs: use only for unit testing of this class
        """
        self._context = context
        self._config_location = test_kwargs.get("config", self._resolve_from_root("config"))
        self._resources_location = test_kwargs.get("resources", self._resolve_from_root("resources"))
        self._dependency_config_location = test_kwargs.get(
            "dependencies", self._resolve_from_root("dependencies/config")
        )
        self._allowed_wf_extensions = {
            "cwl": 1,
            "cwl.json": 2,
            "wdl": 3
        }
        self._suites = {
            "M": "medium.json",
            "L": "large.json",
            "CT": "cwl_tools.json",
            "CW": "cwl_workflows.json",
            "CC": "cwl_cromwell.json"
        }

    def discover(self, suite: str) -> List[RegisteredTest]:
        """
         Method to discover and register the tests.
        Args:
            suite: str. One of several suites supplied as a capital letter code. See self._suites.keys

        Returns: List[RegisteredTest]. List of registered tests
        """
        if suite not in self._suites.keys():
            raise TestDiscoveryError(f"TestDiscovery.discover(): suite {suite} is not registered. Existing suites "
                                     f"are {self._suites}")
        with open(os.path.join(self._config_location, self._suites.get(suite)), "r") as suite_handle:
            suite_config = self._config_linter(json.load(suite_handle))
        registered_tests = [RegisteredTestFactory.register_test(
            src_file=self._find_workflow_source(x[0], x[1]),
            category=x[0],
            test_name=x[1],
            context=self._context
        ) for x in self._flatten_config(suite_config)]
        return registered_tests

    def discover_single_test(self, test_name: str) -> List[RegisteredTest]:
        """
         Method to discover and register a single test. Find the first complete match of the test name.
        Args:
            test_name: str. Test name as present in one of the suite config files

        Returns: List[RegisteredTest]. Subclass of a RegisteredTest, according to the test category. For consistency,
                returned as a list of 1 element
        """
        for suite_file in self._suites.values():
            with open(os.path.join(self._config_location, suite_file), "r") as suite_handle:
                suite_config = self._config_linter(json.load(suite_handle))
            for category, test_collection in suite_config.items():
                if test_name in test_collection:
                    registered_test = RegisteredTestFactory.register_test(
                        src_file=self._find_workflow_source(category, test_name),
                        category=category,
                        test_name=test_name,
                        context=self._context
                    )
                    return [registered_test]
        else:
            raise TestDiscoveryError(
                f"Test name {test_name} was not found among any existing test suites {self._suites}. Add the test to "
                f"the suite. See README `Adding Tests` for instructions"
            )

    def add_tests(self, dir_name: str, extension: str, suite: str, category: str) -> List[str]:
        """
        Adds tests to the test suite.
        Args:
            dir_name: str. Location of the test-related files (workflow, inputs-outputs, etc.)
            extension: str. Test workflow format. Allowed are 'cwl', 'cwl.json', 'wdl'
            suite: str. Test suite name. Usually a team-defined group of tests to be run in each CI/CD step
            category: str. Test category name. Usually a team-defined category that reflects a type of the test

        Returns: List[str]. Names of the added tests
        """
        if suite not in self._suites.keys():
            raise TestDiscoveryError(
                f"TestDiscovery.add_tests(): registering new suite {suite} is not supported. To add a new suite "
                f"consult with the team. Existing suites are {self._suites}"
            )
        added_tests = []
        if extension not in self._allowed_wf_extensions.keys():
            raise TestDiscoveryError(
                f"TestDiscovery.add_tests(): {extension} is not recognized. "
                f"Allowed extensions are {self._allowed_wf_extensions.keys()}"
            )
        for test_file in glob(os.path.join(dir_name, f"*.{extension}")):
            test_name = rm_suffix(os.path.basename(test_file), f".{extension}")
            test_file_with_dependencies = glob(os.path.join(dir_name, f"{test_name}*"))
            self._copy_resource_files(test_file_with_dependencies, category)
            # Next line inefficient I/O, but an attempt to make EVERY test addition atomic
            self._add_name_to_config(test_name, suite, category)
            added_tests.append(test_name)
        logging.warning(f"Added {len(added_tests)} tests to the suite `{suite}` under category `{category}`")
        return added_tests

    def discover_dependencies(self) -> Set[Dependency]:
        """
        Method to discover existing dependencies for every language

        Returns: Set[Dependency]. Set of objects of class Dependency
        """
        dependencies = []
        for dependency_config in os.listdir(self._dependency_config_location):
            if not dependency_config.endswith(".json"):
                continue
            config_path = os.path.join(self._dependency_config_location, dependency_config)
            dependency_factory = DependencyFactory(config_file=config_path)
            one_dependency = dependency_factory.make()
            dependencies.append(one_dependency)
        return set(dependencies)

    @staticmethod
    def _config_linter(config) -> Dict:
        for test_category, tests in config.items():
            if not isinstance(tests, List):
                raise TestDiscoveryError(
                    f"Unexpected config structure for category {test_category}. "
                    f"Expected an array of test names, received {tests}"
                )
            else:
                continue
        return config

    def _resolve_from_root(self, dir_name: str) -> Path:
        return Path(os.path.join(self._context.repo_root_dir, dir_name)).resolve()

    def _add_name_to_config(self, test_name: str, suite: str, category: str) -> None:
        suite_file_path = os.path.join(self._config_location, self._suites.get(suite))
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

    @staticmethod
    def _flatten_config(config: Dict) -> List[Tuple[str, str]]:
        """
        Flattens the config dict.
        Args:
            config: Dict. Config to flatten

        Returns: List[Tuple[str, str]]. Each tuple is Category [0] and test name [1]
        """
        flat_config = []
        for category, test_collection in config.items():
            flattened = [(category, x) for x in test_collection]
            flat_config.extend(flattened)
        return flat_config

    def _find_workflow_source(self, category: str, test_name: str) -> str:
        potential_source_files = []
        for extension in self._allowed_wf_extensions.keys():
            potential_source_files.extend(
                glob(os.path.join(self._resources_location, category, f"{test_name}.{extension}"))
            )
        if len(potential_source_files) == 0:
            raise TestDiscoveryError(
                f"TestDiscovery._find_workflow_source(): could not find any source file for test `{test_name}` in "
                f"category `{category}`. First try to add tests by using `add-tests` CLI command"
            )
        elif len(potential_source_files) > 1:
            raise TestDiscoveryError(
                f"TestDiscovery._find_workflow_source(): for test name `{test_name}` in category `{category}`, "
                f"found source files for multiple workflow languages: {potential_source_files}. "
                f"Delete one of the test files with duplicate name and add again with a changed name."
            )
        else:
            return potential_source_files[0]
