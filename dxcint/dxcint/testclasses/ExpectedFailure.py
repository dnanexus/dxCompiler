from typing import Dict

from dxcint.RegisteredTest import RegisteredTest
from dxcint.Context import Context
from dxcint.Messenger import State


class ExpectedFailure(RegisteredTest):
    def __init__(self, src_file: str, category: str, test_name: str, context: Context):
        super().__init__(src_file, category, test_name, context)

    def _run_executable(self) -> str:
        execution = self._run_executable_inner()
        return execution.describe().get("id")

    def _validate(self) -> Dict:
        self.messenger.wait_for_completion()
        if self.messenger.state == State.FAIL:
            return {
                "passed": True,
                "message": f"Execution of the test {self.name} failed as expected.",
            }
        else:
            return {
                "passed": False,
                "message": f"Execution of the test {self.name} DID NOT fail as expected.",
            }
