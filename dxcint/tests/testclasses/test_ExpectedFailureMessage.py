import os
from dxcint.Messenger import State
from dxcint.testclasses.ExpectedFailureMessage import ExpectedFailureMessage


def test_ExpectedFailureMessage(mock_messenger, fixtures_dir, context_init):
    efm = ExpectedFailureMessage(
        src_file=os.path.join(fixtures_dir, "resources", "mock_category", "efm.wdl"),
        test_name="efm",
        category="mock_category",
        context=context_init,
    )
    efm._messenger = mock_messenger
    efm._messenger.state = State.FAIL
    efm._messenger.describe = {"failureMessage": "expected error"}

    assert efm._validate() == {
        "passed": True,
        "message": "Execution of the test efm failed as expected.",
    }
    efm._messenger.describe = {"failureMessage": "unexpected error"}
    assert efm._validate() == {
        "passed": False,
        "message": "Analysis efm results are invalid.\nExpected: expected error.\nReceived: unexpected error",
    }
    efm._messenger.state = State.FINISHED
    assert efm._validate() == {
        "passed": False,
        "message": "Execution of the test efm DID NOT fail as expected.",
    }
