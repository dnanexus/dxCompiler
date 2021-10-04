import glob
import json
import os
from pathlib import Path
from typing import Optional, Tuple
import urllib.request

from dx_cwl_runner import input_utils, utils
from dx_cwl_runner.dx import Dx


def get_latest_compiler_release():
    with urllib.request.urlopen(
            "https://api.github.com/repos/dnanexus/dxCompiler/releases/latest"
    ) as url:
        return json.loads(url.read().decode())


def download_compiler(release_json) -> Path:
    asset = release_json["assets"][0]
    filename = Path(asset["name"])
    urllib.request.urlretrieve(asset["url"], filename)
    return filename


def pick_latest_compiler_version(jars) -> (Path, str):
    if utils.parse_version:
        versions = dict(
            (utils.parse_version(utils.run_cmd(f"java -jar {jar} version")), jar)
            for jar in jars
        )
        latest_version = list(sorted(versions.keys()))[-1]
        return versions[latest_version], str(latest_version)
    else:
        # pick latest by timestamp
        versions = dict((os.path.getmtime(jar), jar) for jar in jars)
        latest_timestamp = list(sorted(versions.keys()))[-1]
        latest_jar = versions[latest_timestamp]
        latest_version = utils.run_cmd(f"java -jar {latest_jar} version")
        return latest_jar, latest_version


class Compiler:
    def __init__(self, dx: Dx):
        self.dx = dx
        self._init_compiler_jar()
        self._check_tools()

    def _init_compiler_jar(self) -> (Path, str):
        try:
            latest_release = get_latest_compiler_release()
        except Exception as ex:
            self.dx.log.error(f"Error getting latest release from github", ex)
            latest_release = None

        compiler_jar_env: str = os.environ.get("DX_COMPILER_JAR")
        force_update = False
        if compiler_jar_env:
            compiler_jar = Path(compiler_jar_env)
            compiler_version = utils.run_cmd(f"java -jar {compiler_jar} version")
        else:
            jars = glob.glob("dxCompiler*.jar")
            if len(jars) == 1:
                compiler_jar = jars[0]
                compiler_version = utils.run_cmd(f"java -jar {compiler_jar} version")
                force_update = True
            elif len(jars) > 1:
                compiler_jar, compiler_version = pick_latest_compiler_version(jars)
                force_update = True
            elif latest_release:
                compiler_jar = download_compiler(latest_release)
                compiler_version = latest_release["tag"]
                force_update = None
            else:
                raise Exception(
                    "dxCompiler JAR not found in the 'DX_COMPILER_JAR' environment variable, nor "
                    "in the current directory, and the latest version could not be retrieved from GitHub"
                )

        if force_update is not None and latest_release:
            latest_version_str = latest_release["tag_name"]
            if utils.parse_version:
                current_version = utils.parse_version(compiler_version)
                latest_version = utils.parse_version(latest_version_str)
                if current_version < latest_version:
                    if force_update:
                        compiler_jar = download_compiler(latest_release)
                        compiler_version = latest_version_str
                    else:
                        self.dx.log.warning(
                            f"Your version of dxCompiler ({current_version}) is older than the"
                            f"latest release ({latest_version}); you may want to upgrade."
                        )
            elif compiler_version != latest_version_str:
                self.dx.log.warning(
                    f"Your version of dxCompiler is {compiler_version}, but the latest version is "
                    f"{latest_version_str}; you may want to upgrade."
                )

        self.compiler_jar = compiler_jar
        self.compiler_version = compiler_version

    def _check_tools(self):
        pass


class CwlCompiler(Compiler):
    def _check_tools(self):
        utils.check_tool("cwltool", ".+ ([\\d.]+)", "common-workflow-language/cwltool", self.dx.log)
        utils.check_tool("cwl-upgrader", None, "common-workflow-language/cwl-upgrader", self.dx.log)

    def run(self, processfile: str, dx_input_file: str) -> Tuple[Optional[str], Optional[str]]:
        cmd = (
            f"java -jar {self.compiler_jar} compile {processfile} -force "
            f"-folder {self.dx.test_folder} -project {self.dx.current_dx_project} "
            f"-locked -inputs {dx_input_file}"
        )
        if self.dx.log.dryrun:
            self.dx.log.debug(f"compile command: {cmd}")
            return None, None
        else:
            executable = utils.run_cmd(cmd)
            new_dx_input = input_utils.get_new_dx_input(dx_input_file)

            self.dx.log.debug(f"Running {executable} with {new_dx_input} input.")
            job_id = utils.run_cmd(f"dx run {executable} -f {new_dx_input} -y --brief")

            self.dx.log.debug(f"Waiting for {job_id} to finish...")
            utils.run_cmd(f"dx wait {job_id}")
            return job_id, utils.run_cmd(f"dx watch {job_id} --quiet")
