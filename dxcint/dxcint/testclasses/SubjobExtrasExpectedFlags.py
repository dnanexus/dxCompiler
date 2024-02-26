from dxcint.testclasses.ExpectedFlags import ExpectedFlags
from dxcint.mixins.DynamicOnlyMixin import DynamicOnlyMixin
from dxcint.mixins.PinnedInstanceMixin import PinnedInstanceMixin


class SubjobExtrasExpectedFlags(DynamicOnlyMixin, PinnedInstanceMixin, ExpectedFlags):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)