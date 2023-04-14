from dxcint.RegisteredTest import RegisteredTestError
from dxcint.mixins.UnlockedMixin import UnlockedMixin
from dxcint.testclasses.ExpectedOutput import ExpectedOutput
from typing import Dict


class UnlockedExpectedOutput(UnlockedMixin, ExpectedOutput):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _extract_outputs(self) -> Dict:
        desc = self.messenger.describe
        if desc["class"] == "analysis":
            stages = desc["stages"]
            for stage in stages:
                if stage["id"] == "stage-outputs":
                    return stage["execution"]["output"]
            raise RegisteredTestError(
                f"Analysis for test {self.name} does not have stage 'outputs'"
            )
        elif desc["class"] == "job":
            return desc["output"]
        else:
            raise RegisteredTestError(f"Unknown {desc['class']}")
