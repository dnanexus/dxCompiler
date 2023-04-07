import pytest
from dxcint.testclasses.UnlockedExpectedOutput import UnlockedExpectedOutput
from dxcint.Context import ContextEmpty
from dxcint.Messenger import Messenger


@pytest.fixture
def init_LockedExpectedOutput(mocker, context_empty_init):
    mocker.patch.object(
        UnlockedExpectedOutput,
        "job_id",
        return_value="job-XXXX",
        new_callable=mocker.PropertyMock,
    )
    mocker.patch.object(Messenger, "__init__", return_value=None)
    eo = UnlockedExpectedOutput(
        test_name="mock_1", src_file="", category="eo", context=context_empty_init
    )
    yield eo


def test__extract_outputs(init_LockedExpectedOutput, mock_analysis_desc, mocker):
    mocker.patch.object(
        Messenger,
        "describe",
        return_value=mock_analysis_desc,
        new_callable=mocker.PropertyMock,
    )
    outs = init_LockedExpectedOutput._extract_outputs()
    assert outs == {
        "stage-0.done": "Hello World!",
        "stage-common.in_1": "Hello World!",
        "stage-outputs.out": "Hello World!",
    }
