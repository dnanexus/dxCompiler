from dxcint.mixins.ManifestMixin import ManifestMixin
from dxcint.testclasses.AnalysisFinished import AnalysisFinished


class ManifestAnalysisFinished(ManifestMixin, AnalysisFinished):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
