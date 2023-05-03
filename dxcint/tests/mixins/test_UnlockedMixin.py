import os
import subprocess
from dxcint.RegisteredTest import RegisteredTest
from dxcint.mixins.UnlockedMixin import UnlockedMixin


def test_UnlockedMixin(fixtures_dir, context_init, mocker):
    class UnlockedMixer(UnlockedMixin, RegisteredTest):
        pass

    um = UnlockedMixer(
        os.path.join(fixtures_dir, "resources", "mock_category", "mock_1.wdl"),
        "mock_category",
        "mock_1",
        context_init,
    )
    mocker.patch("subprocess.check_output")
    spy = mocker.spy(subprocess, "check_output")
    _ = um.exec_id
    assert "-locked" not in spy.call_args_list[0].args[0]
