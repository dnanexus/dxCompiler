import pytest
from dxcint.testclasses.ExpectedFlags import ExpectedFlags
from dxcint.Messenger import Messenger


@pytest.fixture
def init_ExpectedFlags(mocker, context_empty_init):
    mocker.patch.object(
        ExpectedFlags,
        "job_id",
        return_value="job-XXXX",
        new_callable=mocker.PropertyMock,
    )
    mocker.patch.object(Messenger, "__init__", return_value=None)
    eo = ExpectedFlags(
        test_name="mock_1", src_file="", category="eo", context=context_empty_init
    )
    yield eo


def test__extract_outputs(init_ExpectedFlags, mock_analysis_desc, mocker):
    mocker.patch.object(
        Messenger,
        "describe",
        return_value=mock_analysis_desc,
        new_callable=mocker.PropertyMock,
    )
    outs = init_ExpectedFlags._extract_outputs()
    assert outs == {
        "stage-0.done": "Hello World!",
        "out": "Hello World!",
    }


def test__validate_outputs(
    init_ExpectedFlags, mock_analysis_desc, mock_analysis_results, mocker
):
    mocker.patch.object(
        Messenger,
        "describe",
        return_value=mock_analysis_desc,
        new_callable=mocker.PropertyMock,
    )
    output = init_ExpectedFlags._extract_outputs()
    passed = init_ExpectedFlags._validate_outputs(output, mock_analysis_results)
    assert passed
