from dxcint.mixins.ReorgMixin import ReorgMixin
from dxcint.testclasses.ExpectedOutput import ExpectedOutput


class ReorgExpectedOutput(ReorgMixin, ExpectedOutput):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
