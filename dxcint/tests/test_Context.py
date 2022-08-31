import pytest
import dxpy
from dxpy.api import project_list_folder
from dxcint.Context import ContextError


def test__resolve_project(context_init):
    assert context_init.project_id == list(dxpy.find_projects(name="dxCompiler_playground", level="VIEW"))[0].get("id")


def test__get_version(context_init):
    assert "SNAPSHOT" in context_init._get_version()


def test__create_platform_build_folder(context_init):
    assert "/dxcint_testing/test_" in context_init.platform_build_dir
    list_dir1 = project_list_folder(
        object_id=context_init.project_id,
        input_params={"folder": "/dxcint_testing"}
    )
    last_folder = list_dir1.get("folders")[-1]

    list_dir2 = project_list_folder(
        object_id=context_init.project_id,
        input_params={"folder": last_folder}
    )
    assert f"{last_folder}/applets" in list_dir2.get("folders")


def test_immutability(context_init):
    with pytest.raises(ContextError) as e_info:
        context_init._project_id = "Wrong"
        assert e_info
