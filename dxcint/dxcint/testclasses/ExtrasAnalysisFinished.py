from dxcint.mixins.ExtrasMixin import ExtrasMixin
from dxcint.testclasses.AnalysisFinished import AnalysisFinished


class ExtrasAnalysisFinished(ExtrasMixin, AnalysisFinished):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
