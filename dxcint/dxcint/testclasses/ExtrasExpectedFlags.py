from dxcint.mixins.ExtrasMixin import ExtrasMixin
from dxcint.testclasses.ExpectedFlags import ExpectedFlags


class ExtrasExpectedFlags(ExtrasMixin, ExpectedFlags):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
