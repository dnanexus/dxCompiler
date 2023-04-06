import os
import dxpy
import subprocess
from dxcint.testclasses.ExpectedOutput import ExpectedOutput
from dxcint.RegisteredTest import RegisteredTestError


class ExternExpectedOutput(ExpectedOutput):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._compiler_jar_path = os.path.join(
            self.context.repo_root_dir, f"dxCompiler-{self.context.version}.jar"
        )
        self._applet_folder = os.path.join(
            self._context.platform_build_dir, "extern_applets"
        )
        self._dxni_output_folder = os.path.join(
            self.context.repo_root_dir,
            "dxcint/resources/extern_expected_output/",
        )

    def _compile_executable(self, *args, **kwargs) -> str:
        try:
            with self.context.lock:
                if not os.path.exists(os.path.join(self._dxni_output_folder)):
                    self.context.logger.info("Creating directory for dxni stubs")
                    self._build_native_applets()
                    # works with -folder and with -path
                    self._create_dxni_stub(
                        "-folder",
                        self._applet_folder,
                        os.path.join(self._dxni_output_folder, "dx_extern.wdl"),
                    )
                    self._create_dxni_stub(
                        "-path",
                        os.path.join(self._applet_folder, "/native_concat"),
                        os.path.join(self._dxni_output_folder, "dx_extern_one.wdl"),
                    )
                    self.context.logger.info(f"{self._dxni_output_folder} created")

        except subprocess.CalledProcessError as e:
            self._context.logger.error(
                f"Error compiling {self._dxni_output_folder}\n"
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

    def _create_dxni_stub(self, path_or_folder: str, src: str, dest: str):
        if path_or_folder not in ["-path", "-folder"]:
            raise RegisteredTestError(
                f"path_or_folder argument should be -path or -folder. Provided {path_or_folder}"
            )
        cmdline = [
            "java",
            "-jar",
            self._compiler_jar_path,
            "dxni",
            "-force",
            path_or_folder,
            src,
            "-project",
            self.context.project_id,
            "-language",
            "wdl_v1.0",
            "-output",
            dest,
        ]
        self.context.logger.info(" ".join(cmdline))
        subprocess.check_output(cmdline)
