import os
import pytest
from dxcint.Messenger import State

from dxcint.testclasses.AppExternExpectedOutput import AppExternExpectedOutput
from dxcint.testclasses.ExpectedFailure import ExpectedFailure
from dxcint.testclasses.ExpectedFailureMessage import ExpectedFailureMessage
from dxcint.testclasses.ExpectedOutput import ExpectedOutput
from dxcint.testclasses.ExternExpectedOutput import ExternExpectedOutput
from dxcint.testclasses.ExtrasAnalysisFinished import ExtrasAnalysisFinished
from dxcint.testclasses.LockedExpectedOutput import LockedExpectedOutput
from dxcint.testclasses.ManifestAnalysisFinished import ManifestAnalysisFinished
from dxcint.testclasses.ReorgLockedExpectedOutput import ReorgLockedExpectedOutput
from dxcint.testclasses.StaticDefaultInstanceExpectedOutput import (
    StaticDefaultInstanceExpectedOutput,
)



def test_ExpectedFailureMessage(fixtures_dir, context_init):
    efm = ExpectedFailureMessage(
        src_file=os.path.join(fixtures_dir, "resources", "mock_category", "efm.wdl"),
        test_name="efm",
        category="mock_category",
        context=context_init,
    )
    efm._messenger = FakeMessenger()
    efm._messenger.state = State.FAIL
    efm._messenger.describe_execution = lambda: {"failureMessage": "expected error"}

    assert efm._validate() == {
        "passed": True,
        "message": "Execution of the test efm failed as expected.",
    }
    efm._messenger.describe_execution = lambda: {"failureMessage": "unexpected error"}
    assert efm._validate() == {
        "passed": False,
        "message": "Analysis efm results are invalid.\nExpected: expected error.\nReceived: unexpected error",
    }
    efm._messenger.state = State.FINISHED
    assert efm._validate() == {
        "passed": False,
        "message": "Execution of the test efm DID NOT fail as expected.",
    }

def test_ExternExpectedOutput():
    pytest.fail()

# TODO classes that are not just comibnation of mixins:


def test_AppExternExpectedOutput():
    pytest.fail()





def test_LockedExpectedOutput():
    pytest.fail()
