from dxcint.Context import Context
from dxcint.RegisteredTest import RegisteredTest, RegisteredTestError
from dxcint.testclasses.ReorgLockedExpectedOutput import ReorgLockedExpectedOutput
from dxcint.testclasses.ExpectedFailure import ExpectedFailure
from dxcint.testclasses.ExpectedFailureMessage import ExpectedFailureMessage
from dxcint.testclasses.ExpectedOutput import ExpectedOutput
from dxcint.testclasses.LockedExpectedOutput import LockedExpectedOutput
from dxcint.testclasses.AnalysisFinished import AnalysisFinished
from dxcint.testclasses.ExtrasAnalysisFinished import ExtrasAnalysisFinished
from dxcint.testclasses.ExtrasExpectedOutput import ExtrasExpectedOutput
from dxcint.testclasses.StaticDefaultInstanceExpectedOutput import (
    StaticDefaultInstanceExpectedOutput,
)


class RegisteredTestFactory(object):
    @classmethod
    def register_test(
        cls, src_file: str, category: str, test_name: str, context: Context
    ):
        test_type_switch = {
            "mock_category": RegisteredTest,  # For testing only
            "expected_failure": ExpectedFailure,
            "expected_failure_message": ExpectedFailureMessage,
            "expected_output": ExpectedOutput,
            "analysis_finished": AnalysisFinished,
            "locked_expected_output": LockedExpectedOutput,
            "reorg_locked_expected_output": ReorgLockedExpectedOutput,
            "extras_analysis_finished": ExtrasAnalysisFinished,
            "extras_expected_output": ExtrasExpectedOutput,
            "static_default_instance_expected_output": StaticDefaultInstanceExpectedOutput,
            # ADD NEW CATEGORY HERE
        }
        registered_test = test_type_switch.get(category, None)
        if not registered_test:
            raise RegisteredTestError(
                f"RegisteredTestFactory.register_test(): Category {category} is not recognized. Existing categories "
                f"are {test_type_switch.keys()}. Add new category and new test implementation as a subclass of "
                f"`RegisteredTest`."
            )
        return registered_test(src_file, category, test_name, context)
