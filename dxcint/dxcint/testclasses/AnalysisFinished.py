from typing import Dict

from dxcint.RegisteredTest import RegisteredTest
from dxcint.Context import Context
from dxcint.Messenger import State
from dxcint.utils import async_retry


class AnalysisFinished(RegisteredTest):
    def __init__(self, src_file: str, category: str, test_name: str, context: Context):
        super().__init__(src_file, category, test_name, context)

    @async_retry(max_retries=5)
    def _run_executable(self) -> str:
        execution = self._run_executable_inner()
        return execution.describe().get("id")

    def _validate(self) -> Dict:
        self.messenger.wait_for_completion()
        if self.messenger.state == State.FINISHED:
            return {
                "passed": True,
                "message": f"Execution of the test {self.name} finished successfully as expected.",
            }
        elif self.messenger.state == State.FAIL:
            return {
                "passed": False,
                "message": f"Execution of the test {self.name} FAILED unexpectedly.",
            }
        else:
            return {
                "passed": False,
                "message": f"Execution of the test is in an invalid state: {self.messenger.state}.",
            }
