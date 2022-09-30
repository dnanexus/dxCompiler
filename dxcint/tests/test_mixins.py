import os
import pytest
from dxcint.Messenger import State
from dxcint.RegisteredTest import RegisteredTest

from dxcint.mixins.DefaultInstanceMixin import DefaultInstanceMixin
from dxcint.mixins.ExtrasMixin import ExtrasMixin
from dxcint.mixins.LockedMixin import LockedMixin
from dxcint.mixins.ManifestMixin import ManifestMixin
from dxcint.mixins.ReorgMixin import ReorgMixin
from dxcint.mixins.ResultsTestMixin import ResultsTestMixin
from dxcint.mixins.StaticOnlyMixin import StaticOnlyMixin
from unittest.mock import patch, Mock


@patch("subprocess.check_output")
def test_mixin_composition(context_init):
    """Test that the compile args mixins can be used together"""

    class SuperMixer(
        ExtrasMixin,
        LockedMixin,
        ManifestMixin,
        ReorgMixin,
        StaticOnlyMixin,
        RegisteredTest,
    ):
        pass

    process_mock = Mock()
    attrs = {"strip.return_value": b"workflow-xxxx"}
    process_mock.configure_mock(**attrs)
    hc = SuperMixer("mixer", "mock_category", "mixer", context_init)
    with patch("subprocess.check_output", return_value=process_mock) as mock_subprocess:
        mock = hc._compile_executable()
        compiler_args = mock_subprocess.call_args_list[0].args[0]
        assert "--extras" in compiler_args
        assert "-locked" in compiler_args
        assert "-useManifests" in compiler_args
        assert "-reorg" in compiler_args
        assert "static" in compiler_args
        assert mock == "workflow-xxxx"


# TODO
def test_ExtrasMixin(fixtures_dir, context_init):
    class ExtrasMixer(ExtrasMixin, RegisteredTest):
        pass

    em = ExtrasMixer(
        os.path.join(fixtures_dir, "resources", "mock_category", "mock_1.wdl"),
        "mock_category",
        "mock_1",
        context_init,
    )
    with patch("subprocess.check_output") as mock_subprocess:
        _ = em._compile_executable()
        assert (
            os.path.join(
                fixtures_dir, "resources", "mock_category", "mock_1_extras.json"
            )
            in mock_subprocess.call_args_list[0].args[0]
        )
        assert "--extras" in mock_subprocess.call_args_list[0].args[0]


def test_DefaultInstanceMixin():
    pytest.fail("Not implemented")


def test_ResultsTestMixin():
    pytest.fail("Not implemented")
