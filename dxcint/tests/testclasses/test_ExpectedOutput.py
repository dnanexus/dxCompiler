import pytest
from dxcint.testclasses.ExpectedOutput import ExpectedOutput
from dxcint.Messenger import Messenger


@pytest.fixture
def init_ExpectedOutput(mocker, context_empty_init):
    mocker.patch.object(
        ExpectedOutput,
        "job_id",
        return_value="job-XXXX",
        new_callable=mocker.PropertyMock,
    )
    mocker.patch.object(Messenger, "__init__", return_value=None)
    eo = ExpectedOutput(
        test_name="mock_1", src_file="", category="eo", context=context_empty_init
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
    assert outs == {
        "stage-0.done": "Hello World!",
        "out": "Hello World!",
    }


def test__validate_outputs(
    init_ExpectedOutput, mock_analysis_desc, mock_analysis_results, mocker
):
    mocker.patch.object(
        Messenger,
        "describe",
        return_value=mock_analysis_desc,
        new_callable=mocker.PropertyMock,
    )
    output = init_ExpectedOutput._extract_outputs()
    passed = init_ExpectedOutput._validate_outputs(output, mock_analysis_results)
    assert passed
