import os
import subprocess
from dxcint.RegisteredTest import RegisteredTest
from dxcint.mixins.StaticOnlyMixin import StaticOnlyMixin


def test_StaticOnlyMixin(fixtures_dir, context_init, mocker):
    class StaticOnlyMixer(StaticOnlyMixin, RegisteredTest):
        pass

    som = StaticOnlyMixer(
        os.path.join(fixtures_dir, "resources", "mock_category", "mock_1.wdl"),
        "mock_category",
        "mock_1",
        context_init,
    )
    mocker.patch("subprocess.check_output")
    spy = mocker.spy(subprocess, "check_output")
    _ = som.exec_id
    assert "static" in spy.call_args_list[0].args[0]
    assert "-instanceTypeSelection" in spy.call_args_list[0].args[0]
