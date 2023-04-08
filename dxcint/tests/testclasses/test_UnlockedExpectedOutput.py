import pytest
from dxcint.testclasses.UnlockedExpectedOutput import UnlockedExpectedOutput
from dxcint.Messenger import Messenger


@pytest.fixture
def init_UnlockedExpectedOutput(mocker, context_empty_init):
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


def test__extract_outputs(
    init_UnlockedExpectedOutput, mock_analysis_unlocked_desc, mocker
):
    mocker.patch.object(
        Messenger,
        "describe",
        return_value=mock_analysis_unlocked_desc,
        new_callable=mocker.PropertyMock,
    )
    outs = init_UnlockedExpectedOutput._extract_outputs()
    assert outs == {"out": "Hello World!"}
