from typing import Dict
from dxcint.Messenger import State
from dxcint.Context import Context
from dxcint.testclasses.ExpectedFailure import ExpectedFailure
from dxcint.mixins.ResultsTestMixin import ResultsTestMixin


class ExpectedFailureMessage(ResultsTestMixin, ExpectedFailure):
    def __init__(self, src_file: str, category: str, test_name: str, context: Context):
        super().__init__(src_file, category, test_name, context)
        self._expected_failure_msg = self._get_expected_failure_message()

    def _get_expected_failure_message(self) -> str:
        """Returns the expected failure message from the results."""
        return self.results["error"]

    def _validate(self) -> Dict:
        self.messenger.wait_for_completion()
        if self.messenger.state == State.FAIL:
            failure_msg = self.messenger.execution().describe()["failureMessage"]
            if failure_msg == self._expected_failure_msg:
                return {
                    "passed": True,
                    "message": f"Execution of the test {self.name} failed as expected.",
                }
            else:
                return {
                    "passed": False,
                    "message": f"Analysis {self.name} results are invalid.\nExpected: {self._expected_failure_msg}.\nReceived: {failure_msg}",
                }

        else:
            return {
                "passed": False,
                "message": f"Execution of the test {self.name} DID NOT fail as expected.",
            }
