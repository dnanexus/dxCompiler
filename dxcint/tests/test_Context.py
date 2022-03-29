import pytest
import dxpy
from dxpy.api import project_list_folder
from dxcint.Context import Context, ContextError


@pytest.fixture(scope="session")
def context_init():
    yield Context(project="dxCompiler_playground", folder="/gvaihir/dxcint_testing")


def test__resolve_project(context_init):
    assert context_init.project_id == list(dxpy.find_projects(name="dxCompiler_playground", level="VIEW"))[0].get("id")


def test__get_version(context_init):
    assert "SNAPSHOT" in context_init._get_version()


def test__create_platform_build_folder(context_init):
    assert context_init.platform_build_dir == "/gvaihir/dxcint_testing"
    list_dir = project_list_folder(
        object_id=context_init.project_id,
        input_params={"folder": "/gvaihir/dxcint_testing"}
    )
    assert "/gvaihir/dxcint_testing/applets" in list_dir.get("folders")


def test_immutability(context_init):
    with pytest.raises(ContextError) as e_info:
        context_init._project_id = "Wrong"
        assert e_info
