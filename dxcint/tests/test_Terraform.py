import pytest
import os
from dxcint.Terraform import Terraform
from dxcint.Dependency import DependencyFactory


@pytest.fixture(scope="session")
def terraform_init(dependency_conf_dir, context_init):
    dependency_factory = DependencyFactory(os.path.join(dependency_conf_dir, "mock_dependency.json"))
    package_dependency = dependency_factory.make()
    yield Terraform(languages={"wdl"}, context=context_init, dependencies={package_dependency})


def test___build_compiler(change_to_root_dir, terraform_init):
    compiler_built = terraform_init.build()
    assert compiler_built


def test__wdl_asset(terraform_init):
    built = terraform_init._wdl_asset()
    assert "whatever" in built.keys()


def test__cwl_asset(terraform_init):
    built = terraform_init._cwl_asset()
    assert "whatever" in built.keys()
