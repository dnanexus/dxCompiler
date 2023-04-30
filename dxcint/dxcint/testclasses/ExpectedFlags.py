from typing import Dict
from dxcint.testclasses import ExpectedOutput
from dxcint.mixins.JobCollectorMixin import JobCollectorMixin

from dxpy.api import system_describe_executions


class ExpectedFlags(JobCollectorMixin, ExpectedOutput):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._cache = None

    @property
    def cache(self):
        return self._cache or self._execution_cache()

    def _execution_cache(self) -> Dict:
        exec_ids = self._collect()
        cache = system_describe_executions(exec_ids)
        return {x.get("name"): x for x in cache.get("results", None)}

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
