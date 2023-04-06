import os
import subprocess
from dxcint.RegisteredTest import RegisteredTest
from dxcint.mixins.ReorgMixin import ReorgMixin


def test_ReorgMixin(fixtures_dir, context_init, mocker):
    class ExtrasMixer(ReorgMixin, RegisteredTest):
        pass

    rm = ExtrasMixer(
        os.path.join(fixtures_dir, "resources", "mock_category", "mock_1.wdl"),
        "mock_category",
        "mock_1",
        context_init,
    )
    mocker.patch("subprocess.check_output")
    spy = mocker.spy(subprocess, "check_output")
    _ = rm.exec_id
    assert "-reorg" in spy.call_args_list[0].args[0]
