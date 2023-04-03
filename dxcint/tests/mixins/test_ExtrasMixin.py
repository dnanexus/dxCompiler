import os
import subprocess
from dxcint.RegisteredTest import RegisteredTest
from dxcint.mixins.ExtrasMixin import ExtrasMixin


def test_ExtrasMixin(fixtures_dir, context_init, mocker):
    class ExtrasMixer(ExtrasMixin, RegisteredTest):
        pass

    em = ExtrasMixer(
        os.path.join(fixtures_dir, "resources", "mock_category", "mock_1.wdl"),
        "mock_category",
        "mock_1",
        context_init,
    )
    mocker.patch("subprocess.check_output")
    spy = mocker.spy(subprocess, "check_output")
    _ = em.exec_id
    assert (
        os.path.join(fixtures_dir, "resources", "mock_category", "mock_1_extras.json")
        in spy.call_args_list[0].args[0]
    )
    assert "-extras" in spy.call_args_list[0].args[0]
