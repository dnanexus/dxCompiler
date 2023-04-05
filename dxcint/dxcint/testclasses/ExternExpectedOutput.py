import os
import dxpy
import subprocess
from dxcint.testclasses.ExpectedOutput import ExpectedOutput


class ExternExpectedOutput(ExpectedOutput):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._compiler_jar_path = os.path.join(
            self.context.repo_root_dir, f"dxCompiler-{self.context.version}.jar"
        )
        self._applet_folder = os.path.join(
            self._context.platform_build_dir, "extern_applets"
        )
        self._extern_path = os.path.join(
            self.context.repo_root_dir,
            "dxcint/resources/extern_expected_output/dx_extern.wdl",
        )

    def _compile_executable(self, *args, **kwargs) -> str:
        try:
            with self.context.lock:
                if not os.path.exists(self._extern_path):
                    self.context.logger.info("Creating dx_extern.wdl")
                    self._build_native_applets()
                    self._dxni_create_extern()
                    self.context.logger.info(f"{self._extern_path} created")

        except subprocess.CalledProcessError as e:
            self._context.logger.error(
                f"Error compiling {self._extern_path}\n"
                f"stdout: {e.stdout}\n"
                f"stderr: {e.stderr}"
            )
            self._test_results = {"passed": False, "message": "Error creating extern"}
            raise e
        else:
            return super()._compile_executable(*args, **kwargs)

    def _build_native_applets(self) -> None:
        native_applets = [
            "native_concat",
            "native_diff",
            "native_mk_list",
            "native_sum",
            "native_sum_012",
        ]
        project_handler = dxpy.DXProject(self.context.project_id)
        project_handler.new_folder(self._applet_folder, parents=True)

        # build the native applets, if they do not exist
        for napl in native_applets:
            applet = list(
                dxpy.bindings.search.find_data_objects(
                    classname="applet",
                    name=napl,
                    folder=self._applet_folder,
                    project=self.context.project_id,
                )
            )
            if len(applet) == 0:
                cmdline = [
                    "dx",
                    "build",
                    os.path.join(
                        self.context.repo_root_dir,
                        f"dxcint/dependencies/applets/{napl}",
                    ),
                    "--destination",
                    (self.context.project_id + ":" + self._applet_folder + "/"),
                ]
                self.context.logger.info(" ".join(cmdline))
                subprocess.check_output(cmdline)

    def _dxni_create_extern(self) -> None:
        # build WDL wrapper tasks in test/dx_extern.wdl
        cmdline_common = [
            "java",
            "-jar",
            self._compiler_jar_path,
            "dxni",
            "-force",
            "-folder",
            self._applet_folder,
            "-project",
            self.context.project_id,
        ]

        cmdline_v1 = cmdline_common + [
            "-language",
            "wdl_v1.0",
            "-output",
            self._extern_path,
        ]
        self.context.logger.info(" ".join(cmdline_v1))
        subprocess.check_output(cmdline_v1)
