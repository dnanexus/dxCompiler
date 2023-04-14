import os
import shlex
import random
import dxpy
import json

import subprocess as sp
from typing import Optional, Dict, Union

from dxcint.utils import rm_prefix
from dxcint.Messenger import Messenger
from dxcint.Context import Context, ContextEmpty
from dxcint.constants import DEFAULT_INSTANCE_TYPE, DEFAULT_RUN_KWARGS


class RegisteredTestError(Exception):
    """
    Class to handle RegisteredTest errors
    """


class RegisteredTest(object):
    def __init__(
        self,
        src_file: str,
        category: str,
        test_name: str,
        context: Union[Context, ContextEmpty],
    ):
        self._context = context
        self._src_file = src_file
        self._category = category
        self._test_name = test_name
        self._language = rm_prefix(os.path.basename(src_file), test_name).strip(".")
        self._messenger: Optional[Messenger] = None
        self._exec_id: Optional[str] = None  # id of the executable (workflow or applet)
        self._job_id: Optional[str] = None  # id of the execution (job or analysis)
        self._test_results: Optional[Dict] = None
        self._executable_type_switch = {
            "workflow": dxpy.DXWorkflow,
            "app": dxpy.DXApp,
            "applet": dxpy.DXApplet,
        }
        self._git_revision = sp.check_output(
            ["git", "describe", "--always", "--dirty", "--tags"]
        ).strip()
        # wf inputs supplied as json usually have this suffix. Can be changed in subclasses
        self._raw_inputs_path = os.path.join(
            os.path.dirname(src_file), f"{test_name}_input.json"
        )
        self._compiled_inputs_suffix = "_input.dx.json"
        self._test_inputs = None

    @property
    def context(self) -> Context:
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
    def messenger(self) -> Messenger:
        if not self._messenger:
            self._messenger = self._create_messenger()
        return self._messenger

    @property
    def exec_id(self) -> str:
        if not self._exec_id:
            additional_flags = {"-locked": ""}
            if os.path.exists(self._raw_inputs_path):
                additional_flags = {
                    **additional_flags,
                    **{"-inputs": self._raw_inputs_path},
                }
            self._exec_id = self._compile_executable(
                additional_compiler_flags=additional_flags
            )
        return self._exec_id

    @property
    def job_id(self) -> str:
        if not self._job_id:
            self._job_id = self._run_executable()
        return self._job_id

    @property
    def test_inputs(self) -> Dict:
        if not self._test_inputs:
            self._test_inputs = self._import_inputs()
        return self._test_inputs

    def get_test_result(self) -> bool:
        if not self._test_results:
            self._test_results = self._validate()
        if self._test_results.get("passed"):
            self._context.logger.info(
                f"Test {self._test_name} successfully PASSED."
                f"Additional info from the test: {self._test_results.get('message')}"
            )
            return True
        else:
            self._context.logger.error(
                f"Test {self._test_name} FAILED with message: {self._test_results.get('message')}"
            )
            return False

    def _compile_executable(
        self, additional_compiler_flags: Optional[Dict] = None
    ) -> str:
        """
        Base implementation. For different test classes update this method with new kwargs
        which suite a particular test type. For example, when implementing class ManifestTest(RegisteredTest)
        call this method with `additional_compiler_flags=['-useManifests']` argument in parameter `exec_id`.
        Args:
            additional_compiler_flags: Optional[Dict[str]]. Use this argument to alter the compiler behavior for
            concrete class implementation. Must be provided as dictionary where key is the argument name, value is the
            argument value.

        Returns: str. Compiled workflow ID
        """
        compiler_flags = {
            "-instanceTypeSelection": random.choice(["static", "dynamic"])
        }
        additional_compiler_flags = additional_compiler_flags or {}
        compiler_flags = {
            **compiler_flags,
            **additional_compiler_flags,
        }
        compiler_jar_path = os.path.join(
            self._context.repo_root_dir, f"dxCompiler-{self._context.version}.jar"
        )
        compile_dir = os.path.join(
            self._context.platform_build_dir, "applets", self._test_name
        )
        cmd = (
            f"java -jar {compiler_jar_path} compile {self._src_file} -force -folder {compile_dir} "
            f"-project {self._context.project_id} "
        )
        cmd += " ".join([f"{k} {v}" for k, v in compiler_flags.items()])
        try:
            self._context.logger.info(f"COMPILE COMMAND: {cmd}")
            with self.context.lock:
                workflow_id = sp.check_output(shlex.split(cmd)).strip()
        except sp.CalledProcessError as e:
            self._context.logger.error(
                f"Error compiling {self._src_file}\n"
                f"stdout: {e.stdout}\n"
                f"stderr: {e.stderr}"
            )
            raise e
        return workflow_id.decode("ascii")

    def _run_executable(self, **dx_run_kwargs) -> str:
        """
        This method can be implemented in subclasses with or without additional decorators (e.g. for async_retry)
        Args:
            dx_run_kwargs: Dict[str]. Kwargs for the DxApp(let) or DxWorkflow .run call
        Returns: str. Execution ID (analysis or job)
        """
        execution_desc = self._run_executable_inner(**dx_run_kwargs).describe()
        return execution_desc.get("id")

    def _run_executable_inner(
        self, **dx_run_kwargs
    ) -> Union[dxpy.DXAnalysis, dxpy.DXJob]:
        """
        Override this method if the changes to the workflow execution are needed.

        Returns: Union[DXAnalysis,DXJob]. Execution handler

        """
        if "instance_type" not in dx_run_kwargs:
            dx_run_kwargs.update({"instance_type": DEFAULT_INSTANCE_TYPE})
        dx_run_kwargs.update(DEFAULT_RUN_KWARGS)
        exec_type = self.exec_id.split("-")[0]
        exec_handler = self._executable_type_switch.get(exec_type)(
            project=self._context.project_id, dxid=self.exec_id
        )
        self._context.logger.info(f"Running the process for test {self._test_name}")
        dx_execution = exec_handler.run(
            self.test_inputs,
            project=self._context.project_id,
            folder=os.path.join(self._context.platform_build_dir, "test"),
            name="{} {}".format(self._test_name, self._git_revision),
            **dx_run_kwargs,
        )
        return dx_execution

    def _import_inputs(self) -> Dict:
        input_basename = f"{self._test_name}{self._compiled_inputs_suffix}"
        input_src = os.path.join(os.path.dirname(self._src_file), input_basename)
        if os.path.exists(input_src):
            with open(input_src, "r") as input_src_handle:
                return json.load(input_src_handle)
        else:
            return {}

    def _validate(self) -> Dict:
        """
        Implement this method in subclasses with criteria for test validation (e.g. comparison to output fixtures,
        error messages or other expected behaviors.

        Returns: Dict.
            'passed' - bool. If the test passed or not, according to the implementation
            'message' - str. Auxiliary message specific for the test. Or an error message.
        """
        self.messenger.wait_for_completion()
        return {
            "passed": False,
            "message": "This is an abstract base class. Use a concrete implementation",
        }

    def _create_messenger(self) -> Messenger:
        return Messenger(self.context, self.name, self.job_id)
