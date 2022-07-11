import logging
import os
import shlex
import random
import dxpy

import subprocess as sp
from typing import List, Optional, Dict, Union

from dxcint.utils import rm_prefix
from dxcint.Messenger import Messenger
from dxcint.Context import Context


class RegisteredTestError(Exception):
    """
    Class to handle RegisteredTest errors
    """


class RegisteredTestFactory(object):
    @classmethod
    def register_test(
            cls,
            src_file: str,
            category: str,
            test_name: str,
            context: Context
    ):
        test_type_switch = {
            "mock_category": RegisteredTest,        # For testing only
            # ADD NEW CATEGORY HERE
        }
        registered_test = test_type_switch.get(category, None)
        if not registered_test:
            raise RegisteredTestError(
                f"RegisteredTestFactory.register_test(): Category {category} is not recognized. Existing categories "
                f"are {test_type_switch.keys()}. Add new category and new test implementation as a subclass of "
                f"`RegisteredTest`."
            )
        return registered_test(src_file, category, test_name, context)


class RegisteredTest(object):
    def __init__(
            self,
            src_file: str,
            category: str,
            test_name: str,
            context: Context
    ):
        self._context = context
        self._src_file = src_file
        self._category = category
        self._test_name = test_name
        self._language = rm_prefix(os.path.basename(src_file), test_name).strip(".")
        self._messenger = None
        self._exec_id = None
        self._job_id = None
        self._test_results = None
        self._executable_type_switch = {
            "workflow": dxpy.DXWorkflow,
            "app": dxpy.DXApp,
            "applet": dxpy.DXApplet
        }
        self._git_revision = sp.check_output(
            ["git", "describe", "--always", "--dirty", "--tags"]
        ).strip()

    @property
    def context(self):
        return self._context

    @property
    def category(self) -> str:
        return self._category

    @property
    def name(self) -> str:
        return self._test_name

    @property
    def language(self) -> str:
        return self._language

    @property
    def messenger(self):
        if not self._messenger:
            self._messenger = self._create_messenger()
        return self._messenger

    @property
    def exec_id(self) -> str:
        if not self._exec_id:
            self._exec_id = self._compile_executable()
        return self._exec_id

    @property
    def job_id(self) -> str:
        if not self._job_id:
            self._job_id = self._run_executable()
        return self._job_id

    @property
    def test_result(self) -> bool:
        if not self._test_results:
            self._test_results = self._validate()
        if self._test_results.get("passed"):
            logging.info(
                f"Test {self._test_name} successfully PASSED."
                f"Additional info from the test: {self._test_results.get('message')}"
            )
            return True
        else:
            logging.error(
                f"Test {self._test_name} FAILED with message: {self._test_results.get('message')}"
            )
            return False

    def _compile_executable(self, additional_compiler_flags: Optional[List[str]] = None) -> str:
        """
        Base implementation. For different test classes override `exec_id` property with calling this method with
        arguments which suite a particular test type. For example, when implementing class ManifestTest(RegisteredTest)
        call this method with `additional_compiler_flags=['-useManifests']` argument in parameter `exec_id`.
        Args:
            additional_compiler_flags: Optional[List[str]]. Use this argument to alter the compiler behavior for
            concrete class implementation

        Returns: str. Compiled workflow ID
        """
        compiler_flags = ["-instanceTypeSelection", random.choice(["static", "dynamic"])]
        compiler_flags += additional_compiler_flags or []
        compiler_jar_path = os.path.join(self._context.repo_root_dir, f"dxCompiler-{self._context.version}.jar")
        compile_dir = os.path.join(self._context.platform_build_dir, "applets", self._test_name)
        cmd = f"java -jar {compiler_jar_path} compile {self._src_file} -force -folder {compile_dir} " \
              f"-project {self._context.project_id} "
        cmd += " ".join(compiler_flags)
        try:
            logging.info(f"COMPILE COMMAND: {cmd}")
            workflow_id = sp.check_output(shlex.split(cmd)).strip()
        except sp.CalledProcessError as e:
            logging.error(
                f"Error compiling {self._src_file}\n"
                f"stdout: {e.stdout}\n"
                f"stderr: {e.stderr}"
            )
            raise e
        return workflow_id.decode("ascii")

    def _run_executable(self) -> str:
        """
        This method will be implemented in subclasses with or without additional decorators (e.g. for async_retry)

        Returns: str. Execution ID (analysis or job)
        """
        exec_desc = self._run_executable_inner().describe()
        return exec_desc.get("id")

    def _run_executable_inner(self, exec_input: Optional[Dict] = None) -> Union[dxpy.DXAnalysis, dxpy.DXJob]:
        """
        Override this method if the changes to the workflow execution are needed.
        Args:
            exec_input: Optional[Dict]. Json-ingested inputs for workflow

        Returns: Union[DXAnalysis,DXJob]. Execution handler

        """
        exec_type = self.exec_id.split("-")[0]
        exec_input = exec_input or {}
        exec_handler = self._executable_type_switch.get(exec_type)(
            project=self._context.project_id,
            dxid=self.exec_id
        )
        logging.info(f"Running the process for test {self._test_name}")
        dx_execution = exec_handler.run(
            exec_input,
            project=self._context.project_id,
            folder=os.path.join(self._context.platform_build_dir, "test"),
            name="{} {}".format(self._test_name, self._git_revision)
        )
        return dx_execution

    def _validate(self) -> Dict:
        """
        Implement this method in subclasses with criteria for test validation (e.g. comparison to output fixtures,
        error messages or other expected behaviors.

        Returns: Dict.
            'passed' - bool. If the test passed or not, according to the implementation
            'message' - str. Auxiliary message specific for the test. Or an error message.
        """
        return {
            "passed": False,
            "message": "This is an abstract base class. Use a concrete implementation"
        }


    def _create_messenger(self) -> Messenger:
        raise RegisteredTestError("Messenger is not implemented")
