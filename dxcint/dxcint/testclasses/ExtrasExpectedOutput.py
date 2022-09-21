from dxcint.mixins.ExtrasMixin import ExtrasMixin
from dxcint.testclasses.ExpectedOutput import ExpectedOutput


class ExtrasExpectedOutput(ExtrasMixin, ExpectedOutput):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
