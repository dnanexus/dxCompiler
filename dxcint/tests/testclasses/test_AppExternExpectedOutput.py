import pytest
import subprocess
from dxcint.testclasses.AppExternExpectedOutput import AppExternExpectedOutput
from dxcint.Messenger import Messenger


@pytest.fixture
def init_AppExternExpectedOutput(mocker, context_empty_init):
    mocker.patch.object(
        AppExternExpectedOutput,
        "job_id",
        return_value="job-XXXX",
        new_callable=mocker.PropertyMock,
    )
    mocker.patch.object(Messenger, "__init__", return_value=None)
    eeo = AppExternExpectedOutput(
        test_name="mock_1", src_file="", category="eeo", context=context_empty_init
    )
    yield eeo


def test__create_dxni_app_stub(init_AppExternExpectedOutput, mocker):
    mocker.patch("subprocess.check_output")
    spy = mocker.spy(subprocess, "check_output")
    init_AppExternExpectedOutput._create_dxni_app_stub()
    assert all(
        x in spy.call_args_list[0].args[0]
        for x in [
            "-apps",
            "only",
            "dxcint/resources/app_extern_expected_output/dx_app_extern.wdl",
        ]
    )
