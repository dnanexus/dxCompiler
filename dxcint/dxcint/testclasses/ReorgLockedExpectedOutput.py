from dxcint.mixins.ReorgMixin import ReorgMixin
from dxcint.testclasses.UnlockedExpectedOutput import UnlockedExpectedOutput


class ReorgLockedExpectedOutput(ReorgMixin, UnlockedExpectedOutput):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
