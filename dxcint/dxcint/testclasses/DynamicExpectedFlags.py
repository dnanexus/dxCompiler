from dxcint.testclasses.ExpectedFlags import ExpectedFlags
from dxcint.mixins.DynamicOnlyMixin import DynamicOnlyMixin


class DynamicExpectedFlags(DynamicOnlyMixin, ExpectedFlags):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)