import pytest
from dxcint.testclasses.ExpectedOutput import ExpectedOutput
from dxcint.Context import ContextEmpty
from dxcint.Messenger import Messenger


@pytest.fixture
def init_ExpectedOutput(mocker):
    mocker.patch.object(
        ExpectedOutput,
        "job_id",
        return_value="job-XXXX",
        new_callable=mocker.PropertyMock,
    )
    mocker.patch.object(Messenger, "__init__", return_value=None)
    eo = ExpectedOutput(
        test_name="eotest", src_file="", category="eo", context=ContextEmpty
    )
    yield eo


def test__extract_outputs(init_ExpectedOutput, mock_analysis_desc, mocker):
    mocker.patch.object(
        Messenger,
        "describe",
        return_value=mock_analysis_desc,
        new_callable=mocker.PropertyMock,
    )
    outs = init_ExpectedOutput._extract_outputs()
    assert outs == {"out": "Hello World!"}


def test__compare_result_path(init_ExpectedOutput, mock_analysis_desc, mocker):
    assert False
