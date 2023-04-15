from dxcint.testclasses.UnlockedExpectedOutput import UnlockedExpectedOutput
from dxcint.mixins.ExtrasMixin import ExtrasMixin


class UnlockedExtrasExpectedOutput(ExtrasMixin, UnlockedExpectedOutput):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
