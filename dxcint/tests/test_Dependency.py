import pytest
import os
import shutil
from pathlib import Path
from dxcint.Dependency import Dependency, DependencyFactory, BinaryDependency, PackageDependency, TarballDependency


@pytest.fixture(scope="session")
def link_dest_and_cleanup(dependency_conf_dir, fixtures_dir):
    link_dir = Path(os.path.join(dependency_conf_dir, "links"))

    yield link_dir
    if os.path.exists(link_dir):
        shutil.rmtree(link_dir)


def test_dependency_factory_package(dependency_conf_dir, context_init):
    dependency_factory = DependencyFactory(os.path.join(dependency_conf_dir, "mock_dependency.json"), context_init)
    package_dependency = dependency_factory.make()
    assert isinstance(package_dependency, PackageDependency)
    assert package_dependency.languages == ["wdl", "cwl"]
    assert package_dependency.export_spec() == {
            "name": "mock_dependency",
            "package_manager": "git",
            "version": "1.0.0"
        }


def test_dependency_factory_binary(dependency_conf_dir, link_dest_and_cleanup, context_init):
    dependency_factory = DependencyFactory(os.path.join(dependency_conf_dir, "mock_dependency_src.json"), context_init)
    binary_dependency = dependency_factory.make()
    symlink = os.path.exists(binary_dependency.link(link_dest_and_cleanup))
    assert os.path.exists(symlink)


def test_dependency_factory_tarball(dependency_conf_dir, link_dest_and_cleanup, context_init):
    dependency_factory = DependencyFactory(os.path.join(dependency_conf_dir, "mock_dependency_tar.json"), context_init)
    tarball_dependency = dependency_factory.make()
    assert isinstance(tarball_dependency, TarballDependency)
    assert tarball_dependency.languages == ["wdl", "cwl"]
    symlink = os.path.exists(tarball_dependency.link(link_dest_and_cleanup))
    assert os.path.exists(symlink)
