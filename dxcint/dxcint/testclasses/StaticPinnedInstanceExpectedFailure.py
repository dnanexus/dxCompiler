from dxcint.mixins.StaticOnlyMixin import StaticOnlyMixin
from dxcint.mixins.PinnedInstanceMixin import PinnedInstanceMixin
from dxcint.testclasses.ExpectedFailure import ExpectedFailure


class StaticPinnedInstanceExpectedFailure(
    StaticOnlyMixin, PinnedInstanceMixin, ExpectedFailure
):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
