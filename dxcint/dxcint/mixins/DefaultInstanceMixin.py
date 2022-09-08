from dxcint.RegisteredTest import RegisteredTest
import dxpy
from typing import Union


class DefaultInstanceMixin(RegisteredTest):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _run_executable_inner(self, **kwargs) -> Union[dxpy.DXAnalysis, dxpy.DXJob]:
        return super()._run_executable_inner(instance_type=None)
