import pytest
from dxcint.Terraform import Terraform


@pytest.fixture(scope="session")
def terraform_init():
    yield Terraform()


def build(terraform_init):
    terraform_init.build()
