import os
import subprocess
from dxcint.RegisteredTest import RegisteredTest
from dxcint.mixins.ManifestMixin import ManifestMixin


def test_ManifestMixin(fixtures_dir, context_init, mocker):
    class ExtrasMixer(ManifestMixin, RegisteredTest):
        pass

    mm = ExtrasMixer(
        os.path.join(fixtures_dir, "resources", "mock_category", "mock_1.wdl"),
        "mock_category",
        "mock_1",
        context_init,
    )
    mocker.patch("subprocess.check_output")
    spy = mocker.spy(subprocess, "check_output")
    _ = mm.exec_id
    assert "-useManifests" in spy.call_args_list[0].args[0]
