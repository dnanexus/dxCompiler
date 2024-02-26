import os
import subprocess
from dxcint.RegisteredTest import RegisteredTest
from dxcint.mixins.DynamicOnlyMixin import DynamicOnlyMixin


def test_DynamicOnlyMixin(fixtures_dir, context_init, mocker):
    class DynamicOnlyMixer(DynamicOnlyMixin, RegisteredTest):
        pass

    som = DynamicOnlyMixer(
        os.path.join(fixtures_dir, "resources", "mock_category", "mock_1.wdl"),
        "mock_category",
        "mock_1",
        context_init,
    )
    mocker.patch("subprocess.check_output")
    spy = mocker.spy(subprocess, "check_output")
    _ = som.exec_id
    assert "dynamic" in spy.call_args_list[0].args[0]
    assert "-instanceTypeSelection" in spy.call_args_list[0].args[0]
