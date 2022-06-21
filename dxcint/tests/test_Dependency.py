import pytest
import os
import shutil
from pathlib import Path
from dxcint.Dependency import Dependency, DependencyFactory, BinaryDependency, PackageDependency


@pytest.fixture(scope="session")
def link_dest_and_cleanup(dependency_conf_dir, fixtures_dir):
    link_dir = Path(os.path.join(dependency_conf_dir, "links"))

    yield link_dir
    if os.path.exists(link_dir):
        shutil.rmtree(link_dir)


def test_dependency_factory_package(dependency_conf_dir):
    dependency_factory = DependencyFactory(os.path.join(dependency_conf_dir, "mock_dependency.json"))
    package_dependency = dependency_factory.make()
    assert isinstance(package_dependency, PackageDependency)
    assert package_dependency.languages == ["wdl", "cwl"]
    assert package_dependency.export_spec() == {
            "name": "mock_dependency",
            "package_manager": "git",
            "version": "1.0.0"
        }


def test_dependency_factory_binary(dependency_conf_dir, fixtures_dir, link_dest_and_cleanup):
    dependency_factory = DependencyFactory(os.path.join(dependency_conf_dir, "mock_dependency_src.json"))
    binary_dependency = dependency_factory.make()
    symlink = os.path.exists(binary_dependency.link(link_dest_and_cleanup))
    assert os.path.exists(symlink)
