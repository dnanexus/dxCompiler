import dxpy
from dxcint.mixins.ReorgMixin import ReorgMixin
from dxcint.testclasses.LockedExpectedOutput import LockedExpectedOutput
from typing import Dict


class ReorgLockedExpectedOutput(ReorgMixin, LockedExpectedOutput):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
