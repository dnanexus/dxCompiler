import pytest


def test_build_executable_wdl(build_executable_wdl):
    assert build_executable_wdl.exec_id.startswith("workflow")


def test_init(registered_test_wdl):
    assert registered_test_wdl.category == "mock_category"
    assert registered_test_wdl.name == "mock_1"
    assert registered_test_wdl.language == "wdl"
    assert registered_test_wdl.test_inputs == {"stage-common.in_1": "Hello World!"}


def test_run_wdl(build_executable_wdl):
    job_id = build_executable_wdl.job_id
    assert job_id.startswith("analysis")
