import dxpy
from dxcint.mixins.LockedMixin import LockedMixin
from dxcint.testclasses.ExpectedOutput import ExpectedOutput
from typing import Dict


class LockedExpectedOutput(LockedMixin, ExpectedOutput):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _extract_outputs(self) -> Dict:
        desc = dxpy.describe(self.job_id)
        if desc["class"] in ["analysis", "job"]:
            return desc["output"]
        else:
            raise RuntimeError(f"Unknown {desc['class']}")
