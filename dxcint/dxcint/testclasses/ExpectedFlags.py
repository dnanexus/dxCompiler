from typing import Dict
from dxcint.testclasses.ExpectedOutput import ExpectedOutput


class ExpectedFlags(ExpectedOutput):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _extract_outputs(self) -> Dict:
        return self.messenger.describe
