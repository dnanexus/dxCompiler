import datetime
import pytest
import os
import json
import dxpy
from dxpy.api import container_remove_objects, project_remove_folder
from dxpy.exceptions import DXAPIError
from dxcint.Context import Context
from dxcint.Dependency import DependencyFactory
from dxcint.Terraform import Terraform
from dxcint.RegisteredTest import RegisteredTest
from dxcint.constants import DEFAULT_TEST_PROJECT
from dxcint.Messenger import State


class MockMessenger(object):
    def __init__(self):
        self.wait_for_completion = lambda: None
        self.state = State.NOT_DONE
        self.describe_execution = lambda: None


@pytest.fixture(autouse=True, scope="session")
def prerequisites():
    yield
    available_projects = list(
        dxpy.find_projects(name=DEFAULT_TEST_PROJECT, level="CONTRIBUTE")
    )
    if not available_projects:
        raise ValueError(
            f"Could not find project {DEFAULT_TEST_PROJECT}, you probably need to be logged into the platform"
        )


@pytest.fixture(scope="session")
def root_dir():
    path_link = os.path.join(os.path.dirname(__file__), "../..")
    yield os.path.realpath(path_link)


@pytest.fixture(scope="session")
def fixtures_dir():
    yield os.path.join(os.path.dirname(__file__), "fixtures")


@pytest.fixture(scope="session")
def timestamp():
    yield str(datetime.datetime.now()).replace(" ", "T").replace(":", "-").split(".")[0]


@pytest.fixture(scope="session")
def test_folder(timestamp):
    yield f"/dxcint_testing/test_{timestamp}"


@pytest.fixture(scope="session")
def context_init(root_dir, test_folder):
    context = Context(
        project=DEFAULT_TEST_PROJECT,
        repo_root=root_dir,
        folder=test_folder,
        logger_verbosity="info",
    )
    yield context
    project_remove_folder(
        object_id=context.project_id,
        input_params={"folder": test_folder, "recurse": True, "force": True},
    )


@pytest.fixture(scope="session")
def dependency_conf_dir(fixtures_dir):
    yield os.path.join(fixtures_dir, "dependencies/config/")


@pytest.fixture(scope="session")
def change_to_root_dir(root_dir):
    original_dir = os.path.dirname(__file__)
    os.chdir(root_dir)
    yield
    os.chdir(original_dir)


@pytest.fixture(scope="session")
def terraform_init(dependency_conf_dir, context_init):
    dependency_factory = DependencyFactory(
        config_file=os.path.join(dependency_conf_dir, "mock_dependency_src.json"),
        context=context_init,
    )
    bin_dependency = dependency_factory.make()
    yield Terraform(
        languages={"wdl"}, context=context_init, dependencies={bin_dependency}
    )


@pytest.fixture(scope="session")
def terraform_build(terraform_init):
    yield terraform_init.build()


@pytest.fixture(scope="session")
def registered_test_wdl(fixtures_dir, context_init):
    yield RegisteredTest(
        src_file=os.path.join(fixtures_dir, "resources", "mock_category", "mock_1.wdl"),
        category="mock_category",
        test_name="mock_1",
        context=context_init,
    )


@pytest.fixture(scope="session")
def build_executable_wdl(registered_test_wdl, terraform_build, change_to_root_dir):
    _ = terraform_build
    _ = registered_test_wdl.exec_id
    yield registered_test_wdl
    try:
        container_remove_objects(
            registered_test_wdl.context.project_id,
            {"objects": [registered_test_wdl.exec_id]},
        )
    except DXAPIError as e:
        print(f"Could not find {registered_test_wdl.exec_id}. Exiting")
        print(e)


@pytest.fixture(scope="function")
def mock_messenger():
    yield MockMessenger()


@pytest.fixture(scope="session")
def mock_analysis_desc(fixtures_dir):
    with open(
        os.path.join(fixtures_dir, "mocker", "mock_analysis.json")
    ) as file_handle:
        mock_analysis = json.load(file_handle)
    yield mock_analysis


@pytest.fixture(scope="session")
def mock_analysis_results(fixtures_dir):
    with open(
        os.path.join(fixtures_dir, "resources", "mock_category", "mock_1_results.json")
    ) as file_handle:
        mock_results = json.load(file_handle)
    yield mock_results
