from dxcint.RegisteredTest import RegisteredTest
import dxpy
from typing import Union


class DefaultInstanceMixin(RegisteredTest):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _run_executable_inner(
        self, *args, **kwargs
    ) -> Union[dxpy.DXAnalysis, dxpy.DXJob]:
        kwargs["instance_type"] = None
        return super()._run_executable_inner(*args, **kwargs)
