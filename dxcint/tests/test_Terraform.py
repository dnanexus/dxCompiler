import pytest
from dxcint.Terraform import Terraform


@pytest.fixture(scope="session")
def terraform_init():
    yield Terraform()


def test__wdl_asset(terraform_init):
    built = terraform_init._wdl_asset()
    assert "whatever" in built.keys()


def test__cwl_asset(terraform_init):
    built = terraform_init._cwl_asset()
    assert "whatever" in built.keys()
