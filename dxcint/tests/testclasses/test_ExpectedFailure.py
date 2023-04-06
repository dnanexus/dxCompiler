from dxcint.testclasses.ExpectedFailure import ExpectedFailure
from dxcint.Context import ContextEmpty
from dxcint.Messenger import State


def test_ExpectedFailure(mock_messenger, context_empty_init):
    ef = ExpectedFailure(
        test_name="eftest", src_file="", category="ef", context=context_empty_init
    )
    ef._messenger = mock_messenger
    ef._messenger.state = State.FAIL
    assert ef._validate() == {
        "passed": True,
        "message": "Execution of the test eftest failed as expected.",
    }
    ef._messenger.state = State.FINISHED
    assert ef._validate() == {
        "passed": False,
        "message": "Execution of the test eftest DID NOT fail as expected.",
    }
