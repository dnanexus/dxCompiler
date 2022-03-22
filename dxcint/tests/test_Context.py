import pytest
import dxpy
from dxcint.Context import Context, ContextError


@pytest.fixture(scope="session")
def context_init():
    yield Context(project="dxCompiler_playground")


def test__resolve_project(context_init):
    assert context_init.project_id == list(dxpy.find_projects(name="dxCompiler_playground", level="VIEW"))[0].get("id")


def test_immutability(context_init):
    with pytest.raises(ContextError) as e_info:
        context_init._project_id = "Wrong"
        assert e_info
