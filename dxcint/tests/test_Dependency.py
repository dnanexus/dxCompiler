import pytest
from dxcint.Dependency import Dependency, DependencyFactory, BinaryDependency, PackageDependency


@pytest.fixture(scope="session")
def dependency_init():
    yield Dependency()


def test_(dependency_init):
    assert False
