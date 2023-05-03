import os
from dxcint.RegisteredTest import RegisteredTest
from dxcint.mixins.JobCollectorMixin import JobCollectorMixin


def test_JobCollectorMixin(fixtures_dir, context_init, mocker):
    class JobCollectorMixer(JobCollectorMixin, RegisteredTest):
        pass

    jcm = JobCollectorMixer(
        os.path.join(fixtures_dir, "resources", "mock_category", "mock_1.wdl"),
        "mock_category",
        "mock_1",
        context_init,
    )
    mocker.patch.object(
        JobCollectorMixer,
        "job_id",
        return_value="analysis-GB63Zv00yzZfP82KGkz1F2Qb",
        new_callable=mocker.PropertyMock,
    )
    collected_executions = jcm._collect()
    assert len(collected_executions) == 6
