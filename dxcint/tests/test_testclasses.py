import os
import pytest
from dxcint.Messenger import State

from dxcint.testclasses.AnalysisFinished import AnalysisFinished
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


class FakeMessenger(object):
    def __init__(self):
        self.wait_for_completion = lambda: None
        self.state = State.NOT_DONE
        self.describe_execution = lambda: None


def test_AnalysisFinished():
    af = AnalysisFinished(test_name="aftest", src_file="", category="af", context=None)
    fake_messenger = FakeMessenger()
    fake_messenger.wait_for_completion = lambda: None
    af._messenger = fake_messenger

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
    fake_messenger.state = State.FINISHED
    assert af._validate() == o1
    fake_messenger.state = State.FAIL
    assert af._validate() == o2
    fake_messenger.state = "bad_state"
    assert af._validate() == o3


def test_ExpectedFailure():
    ef = ExpectedFailure(test_name="eftest", src_file="", category="ef", context=None)
    ef._messenger = FakeMessenger()
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


# TODO classes that are not just comibnation of mixins:


def test_AppExternExpectedOutput():
    pytest.fail()


def test_ExternExpectedOutput():
    pytest.fail()


def test_LockedExpectedOutput():
    pytest.fail()
