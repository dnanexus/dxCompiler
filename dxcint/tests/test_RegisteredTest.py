import os.path

import pytest
from dxpy.api import container_remove_objects
from dxpy.exceptions import DXAPIError
from dxcint.RegisteredTest import RegisteredTest


@pytest.fixture(scope="session")
def registered_test_wdl(fixtures_dir):
    yield RegisteredTest(os.path.join(fixtures_dir, "mock_tests", "mock_test.wdl"))


@pytest.fixture(scope="session")
def build_executable_wdl(registered_test_wdl):
    registered_test_wdl.build()
    yield registered_test_wdl
    try:
        container_remove_objects(registered_test_wdl.project_id, {"objects": [registered_test_wdl.exec_id]})
    except DXAPIError as e:
        print(f"Could not find {registered_test_wdl.exec_id}. Exiting")
        print(e)


def test_init(registered_test_wdl):
    assert registered_test_wdl.name == "foo"
    assert registered_test_wdl.kind == "bar"
    assert registered_test_wdl.project_id == "lol"
    assert registered_test_wdl.language == "wdl"
    assert registered_test_wdl.exec_id is None
    assert registered_test_wdl.job_id is None
    assert registered_test_wdl.has_outputs


def test_build_executable_wdl(build_executable_wdl):
    assert build_executable_wdl.exec_id.startswith("workflow")


def test_run_wdl(build_executable_wdl):
    build_executable_wdl.run()
    assert build_executable_wdl.job_id.startswith("analysis")


def test__parse_expected_results(registered_test_wdl):
    outs = registered_test_wdl._parse_expected_results
    assert isinstance(outs, dict)
    assert outs[f"{registered_test_wdl.name}.outs"] == "foo"
    assert isinstance(outs[f"{registered_test_wdl.name}.outs"], str)
