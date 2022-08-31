import os.path

import pytest
from dxpy.api import container_remove_objects
from dxpy.exceptions import DXAPIError
from dxcint.RegisteredTest import RegisteredTest


@pytest.fixture(scope="function")
def registered_test_wdl(fixtures_dir, context_init):
    yield RegisteredTest(
        src_file=os.path.join(fixtures_dir, "resources", "mock_category",  "mock_1.wdl"),
        category="mock_category",
        test_name="mock_1",
        context=context_init
    )


@pytest.fixture(scope="function")
def build_executable_wdl(registered_test_wdl, terraform_init, change_to_root_dir):
    _ = terraform_init.build()
    _ = registered_test_wdl.exec_id
    yield registered_test_wdl
    try:
        container_remove_objects(registered_test_wdl.context.project_id, {"objects": [registered_test_wdl.exec_id]})
    except DXAPIError as e:
        print(f"Could not find {registered_test_wdl.exec_id}. Exiting")
        print(e)


def test_init(registered_test_wdl):
    assert registered_test_wdl.category == "mock_category"
    assert registered_test_wdl.name == "mock_1"
    assert registered_test_wdl.language == "wdl"
    assert registered_test_wdl._test_inputs == {"stage-common.in_1": "Hello World!"}


def test_build_executable_wdl(build_executable_wdl):
    assert build_executable_wdl.exec_id.startswith("workflow")


def test_run_wdl(build_executable_wdl):
    job_id = build_executable_wdl.job_id
    assert job_id.startswith("analysis")

