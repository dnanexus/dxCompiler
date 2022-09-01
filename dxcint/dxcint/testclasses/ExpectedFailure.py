import dxpy

from typing import Dict

from dxcint.RegisteredTest import RegisteredTest
from dxcint.Context import Context


class ExpectedFailure(RegisteredTest):
    def __init__(self, src_file: str, category: str, test_name: str, context: Context):
        super().__init__(src_file, category, test_name, context)

    def _run_executable(self) -> str:
        execution = self._run_executable_inner()
        try:
            execution.wait_on_done()
        except dxpy.DXJobFailureError:
            pass
        return execution.describe().get("id")

    def _validate(self) -> Dict:
        """
        Implement this method in subclasses with criteria for test validation (e.g. comparison to output fixtures,
        error messages or other expected behaviors.

        Returns: Dict.
            'passed' - bool. If the test passed or not, according to the implementation
            'message' - str. Auxiliary message specific for the test. Or an error message.
        """
        execution = self._run_executable_inner()
        try:
            execution.wait_on_done()
            return {
                "passed": False,
                "message": f"Execution of the test {self.name} DID NOT fail as expected."
            }
        except dxpy.DXJobFailureError:
            return {
                "passed": True,
                "message": f"Execution of the test {self.name} failed as expected."
            }
