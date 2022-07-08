import json
import os
import pytest

from dxcint.TestDiscovery import TestDiscovery
from dxcint.RegisteredTest import RegisteredTest
from dxcint.Dependency import Dependency


@pytest.fixture(scope="session")
def add_test_to_existing_category(fixtures_dir, context_init):
    with open(os.path.join(fixtures_dir, "config", "medium.json"), "r") as old_suite_handle:
        old_suite = json.load(old_suite_handle)
    new_tests_dir = os.path.join(fixtures_dir, "new_tests")
    test_discovery = TestDiscovery(
        context=context_init,
        config=os.path.join(fixtures_dir, "config"),
        resources=os.path.join(fixtures_dir, "resources")
    )
    added_tests = test_discovery.add_tests(
        dir_name=new_tests_dir,
        extension="wdl",
        suite="M",
        category="mock_category"
    )
    yield added_tests
    with open(os.path.join(fixtures_dir, "config", "medium.json"), "w") as modified_suite_handle:
        json.dump(old_suite, modified_suite_handle)
    for added_file in os.listdir(new_tests_dir):
        os.remove(os.path.join(fixtures_dir, "resources", "mock_category", added_file))


@pytest.fixture(scope="session")
def add_test_to_new_category(fixtures_dir, context_init):
    with open(os.path.join(fixtures_dir, "config", "medium.json"), "r") as old_suite_handle:
        old_suite = json.load(old_suite_handle)
    new_tests_dir = os.path.join(fixtures_dir, "new_tests")
    test_discovery = TestDiscovery(
        context=context_init,
        config=os.path.join(fixtures_dir, "config"),
        resources=os.path.join(fixtures_dir, "resources")
    )
    added_tests = test_discovery.add_tests(
        dir_name=new_tests_dir,
        extension="wdl",
        suite="M",
        category="new_mock_category"
    )
    yield added_tests
    with open(os.path.join(fixtures_dir, "config", "medium.json"), "w") as modified_suite_handle:
        json.dump(old_suite, modified_suite_handle)
    for added_file in os.listdir(new_tests_dir):
        os.remove(os.path.join(fixtures_dir, "resources", "new_mock_category", added_file))


def test_discover(fixtures_dir, context_init):
    test_discovery = TestDiscovery(
        context=context_init,
        config=os.path.join(fixtures_dir, "config"),
        resources=os.path.join(fixtures_dir, "resources")
    )
    discovered_tests = test_discovery.discover(suite="M")
    assert len(discovered_tests) == 2
    assert all([isinstance(x, RegisteredTest) for x in discovered_tests])
    assert discovered_tests[0].name == "mock_1"
    assert all([x.category == "mock_category" for x in discovered_tests])


def test_discover_single_test(fixtures_dir, context_init):
    test_discovery = TestDiscovery(
        context=context_init,
        config=os.path.join(fixtures_dir, "config"),
        resources=os.path.join(fixtures_dir, "resources")
    )
    discovered_tests = test_discovery.discover_single_test(test_name="mock_1")
    assert len(discovered_tests) == 1
    assert discovered_tests[0].name == "mock_1"
    assert all([x.category == "mock_category" for x in discovered_tests])



def test_add_tests(add_test_to_existing_category, add_test_to_new_category, fixtures_dir):
    assert add_test_to_existing_category[0] == "mock_test"
    with open(os.path.join(fixtures_dir, "config", "medium.json"), "r") as modified_suite_handle:
        modified_suite = json.load(modified_suite_handle)
    assert "mock_test" in modified_suite["mock_category"]
    assert add_test_to_new_category[0] == "mock_test"
    with open(os.path.join(fixtures_dir, "config", "medium.json"), "r") as modified_suite_handle:
        modified_suite = json.load(modified_suite_handle)
    assert "mock_test" in modified_suite["new_mock_category"]
    assert "new_mock_category" in modified_suite.keys()


def test_discover_dependencies(fixtures_dir, context_init):
    test_discovery = TestDiscovery(
        context=context_init,
        dependencies=os.path.join(fixtures_dir, "dependencies/config")
    )
    dependencies = test_discovery.discover_dependencies()
    assert len(dependencies) != 0
    assert all(isinstance(x, Dependency) for x in dependencies)
