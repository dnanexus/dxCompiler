import logging
import os
import shlex
import random
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
            # TODO STOPPED HERE

    def _compile_executable(self, additional_compiler_flags: Optional[List[str]] = None) -> str:
        compiler_flags = ["-instanceTypeSelection", random.choice(["static", "dynamic"])]
        if "manifest" in self._src_file:
            compiler_flags.append("-useManifests")
        compiler_flags += additional_compiler_flags or []
        compiler_jar_path = os.path.join(self._context.repo_root_dir, f"dxCompiler-{self._context.version}.jar")
        compile_dir = os.path.join(self._context.platform_build_dir, "applets", self._test_name)
        cmd = f"java -jar {compiler_jar_path} compile {self._src_file} -force -folder {compile_dir} -project " \
                  f"{self._context.project_id}"
        cmd += " ".join(compiler_flags)
        try:
            logging.info(f"COMPILE COMMAND: {cmd}")
            oid = sp.check_output(shlex.split(cmd)).strip()
        except sp.CalledProcessError as cpe:
            logging.error(
                f"Error compiling {self._src_file}\n"
                f"stdout: {cpe.stdout}\n"
                f"stderr: {cpe.stderr}"
            )
            raise cpe
        return oid.decode("ascii")

    def _create_messenger(self) -> Messenger:
        return Messenger(
            test_name=self.name,
            job_id=self.job_id,
            variant=None  # TODO
        )
