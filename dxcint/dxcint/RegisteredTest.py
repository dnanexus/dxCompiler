import logging
import os
import shlex
import random
import dxpy

import subprocess as sp
from typing import List, Optional

from dxcint.utils import rm_prefix
from dxcint.Messenger import Messenger
from dxcint.Context import Context


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
            context: Context
    ):
        self._context = context
        self._src_file = src_file
        self._category = category
        self._test_name = test_name
        self._language = rm_prefix(os.path.basename(src_file), test_name)
        self._messenger = None
        self._exec_id = None
        self._job_id = None
        self._executable_type_switch = {
            "workflow": dxpy.DXWorkflow,
            "app": dxpy.DXApp,
            "applet": dxpy.DXApplet
        }
        self._git_revision = sp.check_output(
            ["git", "describe", "--always", "--dirty", "--tags"]
        ).strip()

    @property
    def category(self):
        return self._category

    @property
    def name(self):
        return self._test_name

    @property
    def language(self):
        return self._language

    @property
    def messenger(self):
        if not self._messenger:
            self._messenger = self._create_messenger()
        return self._messenger

    @property
    def exec_id(self):
        if not self._exec_id:
            self._exec_id = self._compile_executable()
        return self._exec_id

    @property
    def job_id(self):
        if not self._job_id:
            self._job_id = self._run_executable()
        return self._job_id

    def _compile_executable(self, additional_compiler_flags: Optional[List[str]] = None) -> str:
        """
        Basic implementation. For different test classes override `exec_id` property with calling this method with
        arguments which suite particular test type. For example, when implementing class ManifestTest(RegisteredTest)
        call this method with `additional_compiler_flags=['-useManifests']` argument in parameter `exec_id`.
        Args:
            additional_compiler_flags: Optional[List[str]]. Use this argument to alter the compiler behavior for
            concrete class implementation

        Returns: str.
        """
        compiler_flags = ["-instanceTypeSelection", random.choice(["static", "dynamic"])]
        compiler_flags += additional_compiler_flags or []
        compiler_jar_path = os.path.join(self._context.repo_root_dir, f"dxCompiler-{self._context.version}.jar")
        compile_dir = os.path.join(self._context.platform_build_dir, "applets", self._test_name)
        cmd = f"java -jar {compiler_jar_path} compile {self._src_file} -force -folder {compile_dir} " \
              f"-project {self._context.project_id}"
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

    def _run_executable(self):
        exec_type = self.exec_id.split("-")[0]
        exec_input = {}
        exec_handler = self._executable_type_switch.get(exec_type)(
            project=self._context.project_id,
            dxid=self.exec_id
        )
        dx_execution = exec_handler.run(
            exec_input,
            project=self._context.project_id,
            folder=os.path.join(self._context.platform_build_dir, "test"),
            name="{} {}".format(self._test_name, self._git_revision)
        )
        exec_desc = dx_execution.describe()
        return exec_desc.get("id")

    def _create_messenger(self) -> Messenger:
        raise RegisteredTestError("Messenger is not implemented")
