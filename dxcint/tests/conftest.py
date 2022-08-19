import pytest
import os

from dxcint.Context import Context
from dxcint.Dependency import DependencyFactory
from dxcint.Terraform import Terraform


@pytest.fixture(scope="session")
def root_dir():
    path_link = os.path.join(os.path.dirname(__file__), "../..")
    yield os.path.realpath(path_link)


@pytest.fixture(scope="session")
def fixtures_dir():
    yield os.path.join(os.path.dirname(__file__), "fixtures")


@pytest.fixture(scope="session")
def context_init(root_dir):
    yield Context(
        project="dxCompiler_playground",
        repo_root=root_dir,
        folder="/gvaihir/dxcint_testing",
        logger_verbosity="info"
    )


@pytest.fixture(scope="session")
def dependency_conf_dir(fixtures_dir):
    yield os.path.join(fixtures_dir, "dependencies/config/")


@pytest.fixture(scope="function")
def change_to_root_dir(root_dir):
    original_dir = os.path.dirname(__file__)
    os.chdir(root_dir)
    yield
    os.chdir(original_dir)


@pytest.fixture(scope="function")
def terraform_init(dependency_conf_dir, context_init):
    dependency_factory = DependencyFactory(
        config_file=os.path.join(dependency_conf_dir, "mock_dependency_src.json"),
        context=context_init
    )
    bin_dependency = dependency_factory.make()
    yield Terraform(languages={"wdl"}, context=context_init, dependencies={bin_dependency})
