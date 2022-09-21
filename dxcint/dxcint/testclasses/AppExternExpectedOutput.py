import os
import subprocess
from dxcint.testclasses.ExpectedOutput import ExpectedOutput


class AppExternExpectedOutput(ExpectedOutput):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._compiler_jar_path = os.path.join(
            self.context.repo_root_dir, f"dxCompiler-{self.context.version}.jar"
        )
        self._extern_path = os.path.join(
            self.context.repo_root_dir,
            "dxcint/resources/app_extern_expected_output/dx_app_extern.wdl",
        )

    def _compile_executable(self, *args, **kwargs) -> str:
        try:
            with self.context.lock:
                if not os.path.exists(self._extern_path):
                    self.context.logger.info("Creating dx_app_extern.wdl")
                    self._native_call_app_setup()
                    self.context.logger.info(f"{self._extern_path} created")
        except Exception:
            self.context.logger.error("Error creating dx_app_extern.wdl")
            self._test_results = {False, "Error creating extern"}
        return super()._compile_executable(*args, **kwargs)

    def _native_call_app_setup(self) -> None:
        # build WDL wrapper tasks in test/dx_app_extern.wdl
        cmdline = [
            "java",
            "-jar",
            self._compiler_jar_path,
            "dxni",
            "-apps",
            "only",
            "-force",
            "-language",
            "wdl_v1.0",
            "-output",
            self._extern_path,
        ]
        self.context.logger.info(" ".join(cmdline))
        subprocess.check_output(cmdline)
