from dxcint.mixins.StaticOnlyMixin import StaticOnlyMixin
from dxcint.mixins.PinnedInstanceMixin import PinnedInstanceMixin
from dxcint.testclasses.ExpectedOutput import ExpectedOutput


class StaticPinnedInstanceExpectedOutput(
    StaticOnlyMixin, PinnedInstanceMixin, ExpectedOutput
):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
