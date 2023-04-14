import os
import subprocess
from dxcint.RegisteredTest import RegisteredTest
from dxcint.mixins.ExtrasMixin import ExtrasMixin
from dxcint.mixins.ManifestMixin import ManifestMixin
from dxcint.mixins.PinnedInstanceMixin import PinnedInstanceMixin
from dxcint.mixins.ReorgMixin import ReorgMixin
from dxcint.mixins.ResultsTestMixin import ResultsTestMixin
from dxcint.mixins.StaticOnlyMixin import StaticOnlyMixin


def test_SuperMixer(fixtures_dir, context_init, mocker):
    class SuperMixer(
        ExtrasMixin,
        ManifestMixin,
        PinnedInstanceMixin,
        ReorgMixin,
        ResultsTestMixin,
        StaticOnlyMixin,
    ):
        pass

    sm = SuperMixer(
        os.path.join(fixtures_dir, "resources", "mock_category", "mock_1.wdl"),
        "mock_category",
        "mock_1",
        context_init,
    )
    mocker.patch("subprocess.check_output")
    spy = mocker.spy(subprocess, "check_output")
    _ = sm.exec_id
    assert (
        os.path.join(fixtures_dir, "resources", "mock_category", "mock_1_extras.json")
        in spy.call_args_list[0].args[0]
    )
    flags = [
        "-extras",
        "-useManifests",
        "-locked",
        "-reorg",
        "-instanceTypeSelection",
        "static",
    ]
    assert all(x in spy.call_args_list[0].args[0] for x in flags)
