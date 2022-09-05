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

    async def _validate(self) -> Dict:
        """
        Implement this method in subclasses with criteria for test validation (e.g. comparison to output fixtures,
        error messages or other expected behaviors.

        Returns: Dict.
            'passed' - bool. If the test passed or not, according to the implementation
            'message' - str. Auxiliary message specific for the test. Or an error message.
        """
        self.messenger.wait_for_completion()
        if self.messenger.state == State.FINISHED:
            return {
            "passed": False,
            "message": f"Execution of the test {self.name} DID NOT fail as expected."
            }
        elif self.messenger.state == State.FAIL:
            return {
            "passed": True,
            "message": f"Execution of the test {self.name} failed as expected."
            }
        else:
            return {
            "passed": False,
            "message": f"Execution of the test is in an invalid state: {self.messenger.state}."
            }

