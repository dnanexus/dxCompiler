import dxpy
from typing import Dict
from dxcint.RegisteredTest import RegisteredTestError
from dxcint.testclasses.ExpectedOutput import ExpectedOutput
from dxcint.mixins.JobCollectorMixin import JobCollectorMixin


class ExpectedFlags(JobCollectorMixin, ExpectedOutput):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._cache = None

    @property
    def cache(self):
        return self._cache or self._execution_cache()

    def _execution_cache(self) -> Dict:
        exec_ids = self._collect()
        cache = dxpy.api.system_describe_executions(
            input_params={"executions": exec_ids}
        ).get("results", None)
        if not cache:
            self.context.logger.error(
                f"No executions (jobs/analyses) were found for the test {self.name}"
            )
            raise RegisteredTestError("No executions found")
        return {x["describe"]["name"]: x["describe"] for x in cache}

    def _extract_outputs(self) -> Dict:
        return self.cache

    def _validate_outputs(self, exec_outputs, expected_output) -> bool:
        try:
            for key, expected_val in expected_output.items():
                top_exec_name, stage_name, flag = key.split(".")
                if top_exec_name != self._test_name:
                    self.context.logger.error(
                        "Execution name in output does not match test name"
                    )
                    return False
                if exec_outputs[stage_name][flag] != expected_val:
                    return False
        except Exception:
            return False
        return True
