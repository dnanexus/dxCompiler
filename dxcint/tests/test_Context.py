import pytest
import dxpy
from dxcint.Context import ContextError
from dxcint.constants import DEFAULT_TEST_PROJECT


def test__resolve_project(context_init):
    assert context_init.project_id == list(
        dxpy.find_projects(name=DEFAULT_TEST_PROJECT, level="VIEW")
    )[0].get("id")


def test__get_version(context_init):
    assert "SNAPSHOT" in context_init._get_version()


def test_immutability(context_init):
    with pytest.raises(ContextError) as e_info:
        context_init._project_id = "Wrong"
        assert e_info
