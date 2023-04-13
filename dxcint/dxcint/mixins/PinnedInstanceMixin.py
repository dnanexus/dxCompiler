from dxcint.RegisteredTest import RegisteredTest


class PinnedInstanceMixin(RegisteredTest):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @property
    def job_id(self) -> str:
        if not self._job_id:
            self._job_id = self._run_executable(dx_run_kwargs={"instance_type": None})
        return self._job_id
