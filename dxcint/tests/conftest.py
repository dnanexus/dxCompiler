import pytest
import os

from dxcint.Context import Context


@pytest.fixture(scope="session")
def root_dir():
    path_link = os.path.join(os.path.dirname(__file__), "../..")
    yield os.path.realpath(path_link)


@pytest.fixture(scope="session")
def fixtures_dir():
    yield os.path.join(os.path.dirname(__file__), "fixtures")


@pytest.fixture(scope="session")
def context_init():
    yield Context(project="dxCompiler_playground", folder="/gvaihir/dxcint_testing")


@pytest.fixture(scope="session")
def dependency_conf_dir(fixtures_dir):
    yield os.path.join(fixtures_dir, "dependencies/config/")


@pytest.fixture(scope="function")
def change_to_root_dir(root_dir):
    original_dir = os.path.dirname(__file__)
    os.chdir(root_dir)
    yield
    os.chdir(original_dir)
