import dxpy
from dxcint.mixins.StaticOnlyMixin import StaticOnlyMixin
from dxcint.mixins.DefaultInstanceMixin import DefaultInstanceMixin
from dxcint.testclasses.ExpectedOutput import ExpectedOutput


class StaticDefaultInstanceExpectedOutput(
    StaticOnlyMixin, DefaultInstanceMixin, ExpectedOutput
):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
