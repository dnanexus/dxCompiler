import pytest
import subprocess
from dxcint.testclasses.ExternExpectedOutput import ExternExpectedOutput
from dxcint.Messenger import Messenger


@pytest.fixture
def init_ExternExpectedOutput(mocker, context_empty_init):
    mocker.patch.object(
        ExternExpectedOutput,
        "job_id",
        return_value="job-XXXX",
        new_callable=mocker.PropertyMock,
    )
    mocker.patch.object(Messenger, "__init__", return_value=None)
    eeo = ExternExpectedOutput(
        test_name="mock_1", src_file="", category="eeo", context=context_empty_init
    )
    yield eeo


def test__dxni_look_folder(init_ExternExpectedOutput, mocker):
    mocker.patch("subprocess.check_output")
    spy = mocker.spy(subprocess, "check_output")
    init_ExternExpectedOutput._create_dxni_stub(
        "-folder", "applet_src", "stub_dest.wdl"
    )
    assert all(
        x in spy.call_args_list[0].args[0]
        for x in ["-folder", "applet_src", "stub_dest.wdl"]
    )
