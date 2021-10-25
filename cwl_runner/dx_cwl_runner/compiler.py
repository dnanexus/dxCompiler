from contextlib import contextmanager
import glob
import json
import os
from pathlib import Path
import shutil
import tempfile
from typing import Optional, Tuple
import urllib.request

from dx_cwl_runner import utils
from dx_cwl_runner.dx import Dx
from dx_cwl_runner.utils import Log


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


def get_compiler_jar(download: bool = True, log: utils.Log = Log.get()) -> (Path, str):
    try:
        latest_release = get_latest_compiler_release()
    except Exception as ex:
        log.error(f"Error getting latest release from github", ex)
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
                if download and force_update:
                    compiler_jar = download_compiler(latest_release)
                    compiler_version = latest_version_str
                else:
                    log.warning(
                        f"Your version of dxCompiler ({current_version}) is older than the"
                        f"latest release ({latest_version}); you may want to upgrade."
                    )
        elif compiler_version != latest_version_str:
            log.warning(
                f"Your version of dxCompiler is {compiler_version}, but the latest version is "
                f"{latest_version_str}; you may want to upgrade."
            )

    return compiler_jar, compiler_version


class Compiler:
    def __init__(self, dx: Dx):
        self.dx = dx
        self.compiler_jar, self.compiler_version = get_compiler_jar(log=self.dx.log)
        self._check_tools()
        # use a cache for uploaded files to make sure we upload each file/directory only once
        self._cache = {}

    def _check_tools(self):
        pass

    def upload_file(self, path: str, basedir: str) -> str:
        if path in self._cache:
            return self._cache[path]
        if path.startswith("dx://"):
            return path
        if not os.path.isabs(path):
            path = os.path.join(basedir, path)
        # see if we've already found/uploaded this file
        if not os.path.exists(path):
            raise Exception(f"path does not exist: {path}")
        if not os.path.isfile(path):
            raise Exception(f"not a file {path}")
        return self.dx.find_or_upload_file(path)

    def upload_dir(self, path: str, basedir: str) -> Tuple[str, Optional[str]]:
        if path in self._cache:
            return self._cache[path]
        if path.startswith("dx://"):
            return path, None
        if not os.path.isabs(path):
            path = os.path.join(basedir, path)
        # see if we've already found/uploaded this file
        if not os.path.exists(path):
            raise Exception(f"path does not exist: {path}")
        if not os.path.isdir(path):
            raise Exception(f"not a file {path}")
        dx_uri = self.dx.find_or_upload_dir(path)
        basename = os.path.basename(path)
        return dx_uri, basename

    def write_dx_input(self, new_input: dict, jobfile: str):
        js = json.dumps(new_input, indent=4)
        if self.dx.log.dryrun:
            self.dx.log.log(f"writing input to {jobfile}:\n{js}")
        else:
            with open(jobfile, "w") as dx_input:
                dx_input.write(js)

    @staticmethod
    def get_new_dx_input(input_path: str) -> str:
        basename = os.path.basename(input_path)
        new_input_filename = f"{os.path.splitext(basename)[0]}.dx.json"
        return new_input_filename

    def create_dx_input(
        self, process_file: str, jobfile: str, basedir: str, tmpdir: str
    ) -> Tuple[str, str]:
        """
        Creates the DNAnexus input file. May modify the process file into a form that can be
        compiled by dxCompiler.

        Args:
            process_file: path to the process file
            jobfile: path to the job file
            basedir: inputs base directory
            tmpdir: temporary directory where updated process and jobfiles are to be written

        Returns: (updated_process_file, updated_job_file) - paths to updated process and job files
        """
        pass

    def run(self, processfile: str, dx_input_file: str) -> Optional[dict]:
        cmd = (
            f"java -jar {self.compiler_jar} compile {processfile} -force "
            f"-folder {self.dx.test_folder} -project {self.dx.current_dx_project} "
            f"-locked -inputs {dx_input_file}"
        )
        if self.dx.log.dryrun:
            self.dx.log.debug(f"compile command: {cmd}")
            return None
        else:
            executable = utils.run_cmd(cmd)
            new_dx_input = Compiler.get_new_dx_input(dx_input_file)

            self.dx.log.debug(f"Running {executable} with {new_dx_input} input.")
            exe_id = utils.run_cmd(
                f"dx run {executable} --destination {self.dx.outputs_folder} -f {new_dx_input} -y --brief"
            )

            self.dx.log.debug(f"Waiting for {exe_id} to finish...")
            utils.run_cmd(f"dx wait {exe_id}")
            description = json.loads(utils.run_cmd(f"dx describe {exe_id} --json"))
            execution_log = utils.run_cmd(f"dx watch {exe_id} --quiet")
            if description.get("state") == "done":
                if self.dx.log.verbose and execution_log:
                    self.dx.log.debug(execution_log)
                return description
            else:
                self.dx.log.error(execution_log)
                raise Exception(f"execution of {processfile} failed")

    @staticmethod
    def write_output(file_name: str, results: dict, print_output=True):
        js = json.dumps(results, indent=4)
        with open(file_name, "w") as result_file:
            result_file.write(js)
        if print_output:
            print(js)

    def create_results(self, outputs: dict, process_file: str, outdir: str) -> dict:
        pass

    def create_output(
        self, execution_desc: dict, process_file: str, outdir: str, outfile: str
    ) -> dict:
        """
        Creates the output file from the DNAnexus job/analysis outputs.

        Args:
            execution_desc: The results of describing the execution
            process_file: the process file
            outdir: the directory where output files are to be downloaded
            outfile: the output file where results are to be written

        Returns: the dict of results
        """
        outputs = execution_desc.get("output")
        results = self.create_results(outputs, process_file, outdir)
        Compiler.write_output(outfile, results)
        return results

    @contextmanager
    def _tempdir(self):
        path = tempfile.mkdtemp()
        try:
            yield path
        finally:
            try:
                shutil.rmtree(path)
            except IOError as ex:
                self.dx.log.error(f"Failed to clean up temp dir {path}", ex)

    def run_test(self, processfile, jobfile, basedir, outdir):
        # create a tempdir to write updated CWL and input files
        with self._tempdir() as tmpdir:
            # Convert the jobfile into a dx inputs file; upload any local files/directories.
            # Check whether there are any hard-coded file/directory paths in the CWL; if so,
            # upload them and replace with URIs.
            self.dx.log.debug("Creating DNAnexus inputs...")
            process_file, dx_input_file = self.create_dx_input(processfile, jobfile, basedir, tmpdir)
            # Compile and run the CWL - will throw an exception if the execution fails
            self.dx.log.debug("Compiling and running process file...")
            execution_desc = self.run(process_file, dx_input_file)

        # Convert dx to CWL outputs
        if not self.dx.log.dryrun:
            self.dx.log.debug("Creating outputs...")
            if not os.path.exists(outdir):
                self.dx.log.info(f"Creating output directory {outdir}")
                os.mkdir(outdir)
            elif not os.path.isdir(outdir):
                raise Exception(f"Outdir {outdir} is not a directory!")
            elif not os.access(outdir, os.W_OK):
                raise Exception(f"You need write access to the outdir repository ({outdir})!")
            self.create_output(
                execution_desc,
                process_file,
                outdir,
                f"{os.path.splitext(process_file)[0]}_results.json"
            )

        self.dx.log.debug("Finished successfully!")


class CwlCompiler(Compiler):
    def _check_tools(self):
        utils.check_tool("cwltool", ".+ ([\\d.]+)", "common-workflow-language/cwltool", self.dx.log)
        utils.check_tool("cwl-upgrader", None, "common-workflow-language/cwl-upgrader", self.dx.log)

    def get_modified_input(self, i, basedir: str):
        if type(i) is list:
            return [self.get_modified_input(x, basedir) for x in i]

        # File and Directory are always dicts. There may also be record types.
        # We differentiate based on the value of "class".
        if isinstance(i, dict):
            cls = i.get("class")
            if cls == "File":
                if "location" in i:
                    i["location"] = self.upload_file(i["location"], basedir)
                elif "path" in i:
                    i["location"] = self.upload_file(i.pop("path"), basedir)
                secondary_files = i.get("secondaryFiles")
                if secondary_files is not None:
                    i["secondaryFiles"] = [
                        self.get_modified_input(x, basedir) for x in secondary_files
                    ]
            elif cls == "Directory":
                location = i.get("location", i.pop("path", None))
                if location is not None:
                    location, basename = self.upload_dir(i["location"], basedir)
                    i["location"] = location
                    if "basename" not in i and basename is not None:
                        i["basename"] = basename
                else:
                    listing = i["listing"]
                    if listing is None:
                        raise Exception(
                            f"Directory is missing a 'location', 'path', or 'listing': {i}"
                        )
                    i["listing"] = self.get_modified_input(listing, basedir)
            else:
                i = dict((k, self.get_modified_input(v, basedir)) for k, v in i.items())

        return i

    @staticmethod
    def get_paths_referenced(process_file: str) -> Optional[dict]:
        deps = json.loads(utils.run_cmd(f"cwltool --print-deps {process_file}"))
        return deps.get("secondaryFiles")

    @staticmethod
    def replace_paths(obj, pathmap):
        if isinstance(obj, dict):
            if obj.get("class") in ("File", "Directory"):
                key = obj.get("location", obj.get("path"))
                if key in pathmap:
                    obj["location"] = pathmap[key]
                    return obj
            else:
                return dict((k, CwlCompiler.replace_paths(v, pathmap)) for k, v in obj.items())
        elif isinstance(obj, list):
            return [CwlCompiler.replace_paths(item, pathmap) for item in obj]
        else:
            return obj

    def create_dx_input(
        self, process_file: str, jobfile: str, basedir: str, tmpdir: str
    ) -> Tuple[str, str]:
        # issues:
        #   if file does not exist, exception is thrown and no json is generated, even if some
        #   files were uploaded

        # make sure the file is v1.2

        # pack the file


        # parse the inputs file
        with open(jobfile) as input_file:
            inputs = json.load(input_file)
            process_file = inputs.pop("cwl:tool", process_file)
            dx_inputs = dict(
                (k, self.get_modified_input(v, basedir))
                for k, v in inputs.items()
            )
            dx_input_file = os.path.join(tmpdir, os.path.basename(jobfile))
            self.write_dx_input(dx_inputs, dx_input_file)

        # check if the CWL has any hard-coded paths
        paths_referenced = CwlCompiler.get_paths_referenced(process_file)
        if paths_referenced:
            pathmap = dict(
                (
                    p.get("location", p.get("path")),
                    self.get_modified_input(p, basedir).get("location")
                )
                for p in paths_referenced
                if "location" in p or "path" in p
            )

            with open(process_file, "rt") as inp:
                cwl = json.load(inp)
            dx_cwl = CwlCompiler.replace_paths(cwl, pathmap)
            process_file = os.path.join(tmpdir, os.path.basename(process_file))
            with open(process_file, "wt") as out:
                json.dump(dx_cwl, out)

        return process_file, dx_input_file

    def get_modified_output(self, output: dict, outdir: str) -> dict:
        file_id = next(iter(output.values()))
        output_path = self.dx.get_file_name(outdir, file_id)
        utils.run_cmd(
            f"dx download {next(iter(output.values()))} --output '{output_path}' --no-progress"
        )
        return {
            "class": "File",
            "checksum": str(utils.get_checksum(output_path)),
            "location": output_path,
            "size": os.path.getsize(output_path),
        }

    def create_results(self, outputs: dict, process_file: str, outdir: str) -> dict:
        results = {}
        for output in outputs:
            if type(outputs[output]) is list:
                results[f"{process_file}.{output}"] = []
                for idx, value in enumerate(outputs[output]):
                    results[f"{process_file}.{output}"][idx] = self.get_modified_output(
                        outputs[outputs][idx], outdir
                    )
            if type(outputs[output]) is dict:
                results[f"{process_file}.{output}"] = self.get_modified_output(
                    outputs[output], outdir
                )
            else:
                results[f"{process_file}.{output}"] = outputs[output]
                continue
        return results
