import dxpy
from dxcint.mixins.ExtrasMixin import ExtrasMixin
from dxcint.testclasses.ExpectedOutput import ExpectedOutput
from typing import Dict


class ExtrasExpectedOutput(ExtrasMixin, ExpectedOutput):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
