import os
from dxcint.RegisteredTest import RegisteredTest
from dxcint.mixins.PinnedInstanceMixin import PinnedInstanceMixin


def test_PinnedInstanceMixin(fixtures_dir, context_init, mocker):
    class PinnedInstnceMixer(PinnedInstanceMixin, RegisteredTest):
        pass

    pim = PinnedInstnceMixer(
        os.path.join(fixtures_dir, "resources", "mock_category", "mock_1.wdl"),
        "mock_category",
        "mock_1",
        context_init,
    )
    mocker.patch.object(RegisteredTest, "_run_executable")
    spy = mocker.spy(PinnedInstanceMixin, "_run_executable")
    _ = pim.job_id
    assert "instance_type" in spy.call_args_list[0].kwargs
    assert spy.call_args_list[0].kwargs["instance_type"] is None
