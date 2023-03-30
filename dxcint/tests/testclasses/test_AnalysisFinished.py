from dxcint.testclasses.AnalysisFinished import AnalysisFinished
from dxcint.Context import ContextEmpty
from dxcint.Messenger import State


def test_AnalysisFinished(mock_messenger):
    af = AnalysisFinished(
        test_name="aftest", src_file="", category="af", context=ContextEmpty
    )
    mock_messenger.wait_for_completion = lambda: None
    af._messenger = mock_messenger

    o1 = {
        "passed": True,
        "message": "Execution of the test aftest finished successfully as expected.",
    }
    o2 = {
        "passed": False,
        "message": "Execution of the test aftest FAILED unexpectedly.",
    }
    o3 = {
        "passed": False,
        "message": "Execution of the test is in an invalid state: bad_state.",
    }
    mock_messenger.state = State.FINISHED
    assert af._validate() == o1
    mock_messenger.state = State.FAIL
    assert af._validate() == o2
    mock_messenger.state = "bad_state"
    assert af._validate() == o3
