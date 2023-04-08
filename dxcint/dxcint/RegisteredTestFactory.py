from dxcint.Context import Context
from dxcint.RegisteredTest import RegisteredTest, RegisteredTestError
from dxcint.testclasses.ReorgExpectedOutput import ReorgExpectedOutput
from dxcint.testclasses.ExpectedFailure import ExpectedFailure
from dxcint.testclasses.ExpectedFailureMessage import ExpectedFailureMessage
from dxcint.testclasses.ExpectedOutput import ExpectedOutput
from dxcint.testclasses.UnlockedExpectedOutput import UnlockedExpectedOutput
from dxcint.testclasses.AnalysisFinished import AnalysisFinished
from dxcint.testclasses.ExtrasAnalysisFinished import ExtrasAnalysisFinished
from dxcint.testclasses.ExtrasExpectedOutput import ExtrasExpectedOutput
from dxcint.testclasses.StaticPinnedInstanceExpectedOutput import (
    StaticPinnedInstanceExpectedOutput,
)
from dxcint.testclasses.ManifestAnalysisFinished import ManifestAnalysisFinished
from dxcint.testclasses.ExternExpectedOutput import ExternExpectedOutput
from dxcint.testclasses.AppExternExpectedOutput import AppExternExpectedOutput


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
            "unlocked_expected_output": UnlockedExpectedOutput,
            "reorg_expected_output": ReorgExpectedOutput,
            "extras_analysis_finished": ExtrasAnalysisFinished,
            "extras_expected_output": ExtrasExpectedOutput,
            "static_pinned_instance_expected_output": StaticPinnedInstanceExpectedOutput,
            "manifest_analysis_finished": ManifestAnalysisFinished,
            "extern_expected_output": ExternExpectedOutput,
            "app_extern_expected_output": AppExternExpectedOutput,
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
