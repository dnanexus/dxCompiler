import dxpy
from typing import Dict
from dxcint.RegisteredTest import RegisteredTestError
from dxcint.testclasses.ExpectedOutput import ExpectedOutput
from dxcint.mixins.JobCollectorMixin import JobCollectorMixin


class ExpectedFlags(JobCollectorMixin, ExpectedOutput):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._cache = {}

    @property
    def cache(self):
        return self._cache or self._execution_cache()

    def _execution_cache(self) -> Dict:
        relevant_flags = self._collect_relevant_flags()
        exec_ids = self._collect()
        describe_payload = [{"id": exec_id, "describe": {"fields": relevant_flags}} for exec_id in exec_ids]
        cache = dxpy.api.system_describe_executions(
            input_params={"executions": describe_payload}
        ).get("results", None)
        if not cache:
            self.context.logger.error(
                f"No executions (jobs/analyses) were found for the test {self.name}"
            )
            raise RegisteredTestError("No executions found")
        for executable in cache:
            # differentiate sub-jobs
            body_suffix = "body" if executable["describe"]["function"] == "body" else ""
            exec_name = executable["describe"]["executableName"]
            full_exec_name = ":".join([
                x for x in [exec_name, body_suffix] if x
            ])
            self._cache.update(**{full_exec_name: executable["describe"]})
        return self._cache

    def _extract_outputs(self) -> Dict:
        return self.cache

    def _validate_outputs(self, exec_outputs, expected_output) -> bool:
        try:
            for key, expected_val in expected_output.items():
                top_exec_name, stage_name, flag = key.split(".")
                if top_exec_name != self.name:
                    self.context.logger.error(
                        "Execution name in output does not match test name"
                    )
                    return False
                if exec_outputs[stage_name][flag] != expected_val:
                    return False
        except Exception:
            return False
        return True

    def _collect_relevant_flags(self) -> Dict:
        """
        Collects all flags in the expected results. For each executable, all collected fags will be requested in the API call
        :return: Dict[String, Bool]
        """
        bag_of_flags = {
            "executableName": True,
            "function": True
        }
        for key in self._results.keys():
            _, _, flag = key.split(".")
            bag_of_flags.update(**{flag: True})
        return bag_of_flags
