import os
import json
import dxpy
from dxpy.api import project_list_folder


def test__build_compiler(change_to_root_dir, terraform_init):
    compiler_built = terraform_init._build_compiler()
    assert compiler_built


def test__create_platform_build_folder(terraform_init, test_folder):
    assert terraform_init._create_build_subdirs()
    list_dir = project_list_folder(
        object_id=terraform_init.context.project_id,
        input_params={"folder": test_folder},
    )
    assert f"{test_folder}/applets" in list_dir.get("folders")


def test__generate_config_file(change_to_root_dir, terraform_init):
    conf = terraform_init._generate_config_file()
    assert 'region = "aws:us-east-1"' in conf
    assert terraform_init._context.platform_build_dir in conf


def test__create_local_asset_dir(change_to_root_dir, terraform_init):
    _ = terraform_init._create_local_asset_dir("WDL")
    assert os.path.exists(
        terraform_init._local_created_dirs.get("WDL").get("resources")
    )


def test__create_asset_spec(change_to_root_dir, terraform_init):
    _ = terraform_init._create_local_asset_dir("WDL")
    spec = terraform_init._create_asset_spec("WDL")
    spec_location = os.path.join(
        terraform_init._context.repo_root_dir, "applet_resources", "WDL", "dxasset.json"
    )
    with open(spec_location, "r") as spec_handle:
        spec_import = json.load(spec_handle)
    assert spec == spec_import
    assert None not in spec_import.get("execDepends")


def test_build(terraform_build, change_to_root_dir, context_init):
    built_ids = terraform_build
    assert len(built_ids) == 1
    descr = dxpy.describe(built_ids[0])
    assert descr.get("name", None) == "dxWDLrt"
    assert descr.get("folder", None) == context_init.platform_build_dir
