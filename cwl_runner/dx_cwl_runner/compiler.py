from contextlib import contextmanager
import glob
import json
import os
from pathlib import Path
import re
import shutil
import tempfile
from typing import Optional, Tuple
import urllib.request
from urllib.parse import quote as urlquote

from dx_cwl_runner import utils
from dx_cwl_runner.dx import Dx


_latest_release = None


def get_latest_compiler_release():
    global _latest_release
    if _latest_release is None:
        with urllib.request.urlopen(
                "https://api.github.com/repos/dnanexus/dxCompiler/releases/latest"
        ) as url:
            _latest_release = json.loads(url.read().decode())
    return _latest_release


def download_compiler(release_json) -> Path:
    asset = release_json["assets"][0]
    filename = Path(asset["name"])
    urllib.request.urlretrieve(asset["browser_download_url"], filename)
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


def get_compiler_jar(log, download: bool = True) -> (Path, str):
    try:
        latest_release = get_latest_compiler_release()
    except:
        log.error(f"Error getting latest release from github", exc_info=True)
        latest_release = None

    compiler_jar_env: str = os.environ.get("DX_COMPILER_JAR")
    force_update = False
    if compiler_jar_env:
        compiler_jar = Path(compiler_jar_env)
        compiler_version = utils.run_cmd(f"java -jar {compiler_jar} version", log)
    else:
        jars = glob.glob("dxCompiler*.jar")
        if len(jars) == 1:
            compiler_jar = jars[0]
            compiler_version = utils.run_cmd(f"java -jar {compiler_jar} version", log)
            force_update = True
        elif len(jars) > 1:
            compiler_jar, compiler_version = pick_latest_compiler_version(jars)
            force_update = True
        elif latest_release:
            if download:
                compiler_jar = download_compiler(latest_release)
            else:
                compiler_jar = None
            compiler_version = latest_release["tag_name"]
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


_folder_uri_re = re.compile("dx://(.+):(.+)")


class Compiler:
    def __init__(self, dx: Dx):
        self.dx = dx
        self.compiler_jar, self.compiler_version = get_compiler_jar(self.dx.log)
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

    def write_modified_input(self, new_input: dict, jobfile: str):
        js = json.dumps(new_input, indent=4)
        if self.dx.log.dryrun:
            self.dx.log.info(f"writing input to {jobfile}:\n{js}")
        else:
            with open(jobfile, "w") as dx_input:
                dx_input.write(js)

    def create_compiler_input(
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

    def run(self, process_file: str, input_file: str) -> Optional[dict]:
        cmd = (
            f"java -jar {self.compiler_jar} compile {process_file} -force "
            f"-folder {self.dx.test_folder} -project {self.dx.current_dx_project} "
            f"-locked -inputs {input_file}"
        )
        if self.dx.log.dryrun:
            self.dx.log.info(f"Compile command: {cmd}")
            return None
        else:
            executable = utils.run_cmd(cmd, self.dx.log)
            self.dx.log.debug(f"Compiled {process_file} to {executable}")
            new_dx_input_file = os.path.join(
                os.path.dirname(input_file),
                f"{os.path.splitext(os.path.basename(input_file))[0]}.dx.json"
            )
            self.dx.log.debug(f"Running {executable} with {new_dx_input_file} input.")
            with open(new_dx_input_file, "r") as inp:
                new_dx_input = json.load(inp)
            description, execution_log = self.dx.run_executable(executable, new_dx_input)
            if description.get("state") == "done":
                if self.dx.log.verbose and execution_log:
                    self.dx.log.debug(execution_log)
                return description
            else:
                self.dx.log.error(execution_log)
                raise Exception(f"execution of {process_file} failed")

    def create_results(self, outputs: dict, outdir: str) -> dict:
        pass

    def create_output(
        self, execution_desc: dict, process_file: str, outdir: str
    ) -> dict:
        """
        Creates the output file from the DNAnexus job/analysis outputs.

        Args:
            execution_desc: The results of describing the execution
            process_file: the process file
            outdir: the directory where output files are to be downloaded

        Returns: the dict of results
        """
        outputs = execution_desc.get("output")
        results = self.create_results(outputs, outdir)
        outfile = os.path.join(
            outdir, f"{os.path.splitext(os.path.basename(process_file))[0]}_results.json"
        )
        with open(outfile, "w") as result_file:
            json.dump(results, result_file, indent=4)
        return results

    @contextmanager
    def _tempdir(self):
        path = tempfile.mkdtemp()
        try:
            yield path
            try:
                shutil.rmtree(path)
            except IOError:
                self.dx.log.error(f"Failed to delete up temp dir {path}", exc_info=True)
        except:
            # don't delete the tmpdir if there is an exception, for debugging purposes
            self.dx.log.warning(
                f"There was an exception - not deleting temporary files in {path}"
            )
            raise

    def run_test(self, process_file, job_file, basedir, outdir):
        def resolve(path, name):
            if not os.path.exists(path) and not os.path.isabs(path) and basedir:
                path = os.path.join(basedir, path)
                if not os.path.exists(path):
                    raise Exception(f"{name} {path} does not exist")
            return path

        process_file = resolve(process_file, "process file")
        job_file = resolve(job_file, "job file")

        # create a tempdir to write updated CWL and input files
        with self._tempdir() as tmpdir:
            # Convert the jobfile into a dx inputs file; upload any local files/directories.
            # Check whether there are any hard-coded file/directory paths in the CWL; if so,
            # upload them and replace with URIs.
            self.dx.log.debug("Creating DNAnexus inputs...")
            modified_process_file, modified_job_file = self.create_compiler_input(
                process_file, job_file, basedir, tmpdir
            )
            # Compile and run the CWL - will throw an exception if the execution fails
            self.dx.log.debug("Compiling and running process file...")
            execution_desc = self.run(modified_process_file, modified_job_file)
        # Convert dx to CWL outputs
        if not self.dx.log.dryrun:
            self.dx.log.debug("Creating outputs...")
            if not os.path.exists(outdir):
                self.dx.log.debug(f"Creating output directory {outdir}")
                os.mkdir(outdir)
            elif not os.path.isdir(outdir):
                raise Exception(f"Outdir {outdir} is not a directory!")
            elif not os.access(outdir, os.W_OK):
                raise Exception(f"You need write access to the outdir repository ({outdir})!")
            self.create_output(
                execution_desc,
                modified_process_file,
                outdir
            )

        self.dx.log.debug("Finished successfully!")


def get_paths_referenced(process_file: str, basedir: Optional[str] = None) -> Optional[dict]:
    deps = json.loads(utils.run_cmd(
        f"cwltool --print-deps --relative-deps cwd {process_file}", cwd=basedir
    ))
    return deps.get("secondaryFiles")


def as_file_uri(x: str, basedir: Optional[str] = None) -> str:
    if x.startswith("file:"):
        return x
    else:
        if not os.path.isabs(x):
            if basedir:
                x = os.path.join(basedir, x)
            else:
                x = os.path.abspath(x)
        return f"file://{urlquote(os.path.realpath(x))}"


def replace_paths(obj, pathmap):
    if isinstance(obj, dict):
        if obj.get("class") in ("File", "Directory"):
            key = obj.get("location", obj.get("path"))
            if key in pathmap:
                obj["location"] = pathmap[key]
            return obj
        else:
            return dict((k, replace_paths(v, pathmap)) for k, v in obj.items())
    elif isinstance(obj, list):
        return [replace_paths(item, pathmap) for item in obj]
    else:
        return obj


def get_main_process(cwl: dict) -> dict:
    if cwl.get("id") == "#main":
        return cwl
    elif "$graph" in cwl:
        for process in cwl["$graph"]:
            if process.get("id") == "#main":
                return cwl
    raise Exception(f"cannot file main process in {cwl}")


def get_schema_defs(proc: dict) -> dict:
    schema_defs = {}
    for req in (proc.get("requirements", []) + proc.get("hints", [])):
        if req["class"] == "SchemaDefRequirement":
            for t in req["types"]:
                schema_defs[t["name"]] = t
    return schema_defs


def dir_to_listing(dir) -> list:
    def entry_to_obj(entry: os.DirEntry):
        if entry.is_file():
            return {
                "class": "File",
                "location": entry.path,
                "basename": entry.name,
                "size": entry.stat().st_size,
                "checksum": str(utils.get_checksum(entry.path)),
            }
        elif entry.is_dir():
            return {
                "class": "Directory",
                "location": entry.path,
                "basename": entry.name,
                "listing": dir_to_listing(entry.path),
            }
        else:
            raise Exception(f"not a file or directory: {entry.path}")

    return [entry_to_obj(entry) for entry in os.scandir(dir)]


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

    # def expand_patterns(self, i: dict, value, schema_defs: dict):
    #     def get_schema_type(t):
    #         if isinstance(t, str):
    #             schema_defs.get(t)
    #         elif isinstance(t, dict) and t.get("type") == "array":
    #             get_schema_type(t["items"])
    #         elif isinstance(t, list):
    #             schema_types = list(filter(None, [get_schema_type(x) for x in t]))
    #             if len(schema_types) == 1:
    #                 return schema_types[0]
    #             elif len(schema_types) > 1:
    #                 raise Exception("more than one schema type")
    #
    #     if (
    #         isinstance(value, dict)
    #         and value.get("class") in {"File", "Directory"}
    #         and "secondaryFiles" not in value
    #     ):
    #         # we only care about file inputs that don't already have secondaryFiles defined
    #         if "secondaryFiles" in i:
    #             for sf in i["secondaryFiles"]:
    #                 pass # TODO
    #         else:
    #             t = get_schema_type(i.get("type"))
    #             if t:
    #                 pass # TODO
    #     else:
    #         return value

    def create_compiler_input(
        self, process_file: str, job_file: str, basedir: str, tmpdir: str
    ) -> Tuple[str, str]:
        # issues:
        #   if file does not exist, exception is thrown and no json is generated, even if some
        #   files were uploaded

        # parse the inputs file
        with open(job_file) as input_file:
            inputs = json.load(input_file)
            if "cwl:tool" in inputs:
                tool = inputs.pop("cwl:tool")
                if process_file is None:
                    process_file = tool
            modified_inputs = dict(
                (k, self.get_modified_input(v, basedir))
                for k, v in inputs.items()
            )

        # try to pack the file
        packed_file = os.path.join(tmpdir, os.path.basename(process_file))
        if not packed_file.endswith(".json"):
            packed_file = f"{packed_file}.json"
        utils.run_cmd(
            f"cwltool --pack {process_file} > {packed_file}", self.dx.log, verbose=False
        )

        # make sure the file is v1.2
        try:
            utils.run_cmd(f"cwl-upgrader {packed_file}", self.dx.log, verbose=False)
        except:
            # it's already 1.2
            pass

        # read the CWL
        with open(packed_file, "rt") as inp:
            cwl = json.load(inp)

        # expand any patterns in input secondaryFiles
        # TODO: we may not need this if cwltool can implement an option to print this information
        #  https://github.com/common-workflow-language/cwltool/issues/1556
        # main_process = get_main_process(cwl)
        # schema_defs = get_schema_defs(cwl)
        # for i in main_process.get("inputs", []):
        #     assert i["id"].startswith("#main/")
        #     name = i["id"][6:]
        #     if name in modified_inputs:
        #         modified_inputs[name] = self.expand_patterns(i, modified_inputs[name], schema_defs)

        # write out the modified jobfile
        modified_jobfile = os.path.join(tmpdir, os.path.basename(job_file))
        self.write_modified_input(modified_inputs, modified_jobfile)

        # if the CWL has any hard-coded local paths, replace them with platform paths
        paths_referenced = get_paths_referenced(packed_file, basedir)
        if paths_referenced:
            pathmap = dict(
                (
                    as_file_uri(p.get("location", p.get("path")), basedir),
                    self.get_modified_input(p, basedir).get("location")
                )
                for p in paths_referenced
                if "location" in p or "path" in p
            )
            # update the packed file in-place
            with open(packed_file, "wt") as out:
                json.dump(replace_paths(cwl, pathmap), out)

        return packed_file, modified_jobfile

    def get_modified_output(self, output, outdir: str):
        if isinstance(output, list):
            return [self.get_modified_output(item, outdir) for item in output]

        if isinstance(output, dict):
            if "___" in output and len(output) == 1:
                output = output["___"]

            cls = output.get("type")
            if cls == "File":
                link = output["uri"]["$dnanexus_link"]
                if isinstance(link, str):
                    file_id = link
                else:
                    file_id = link["id"]
                if "basename" in output:
                    filename = output["basename"]
                else:
                    filename = self.dx.get_file_name(file_id)
                output_path = utils.get_unique_path(outdir, filename)
                self.dx.download_file(file_id, output_path)
                return {
                    "class": "File",
                    "location": output_path,
                    "basename": os.path.basename(output_path),
                    "checksum": output.get("checksum") or str(utils.get_checksum(output_path)),
                    "size": output.get("size") or os.path.getsize(output_path),
                }
            elif cls == "Folder":
                project, folder = _folder_uri_re.match(output["uri"]).groups()
                basename = output.get("basename") or os.path.basename(folder)
                folder_outdir = utils.get_unique_path(outdir, basename)
                if "listing" in output:
                    listing = [
                        self.get_modified_output(x, folder_outdir)
                        for x in output["listing"]
                    ]
                else:
                    self.dx.download_folder(folder, folder_outdir, project)
                    listing = dir_to_listing(folder_outdir)
                return {
                    "class": "Directory",
                    "location": folder_outdir,
                    "basename": basename,
                    "listing": listing
                }
            elif cls == "Listing":
                basename = output["basename"]
                listing_outdir = utils.get_unique_path(outdir, basename)
                return {
                    "class": "Directory",
                    "location": listing_outdir,
                    "basename": basename,
                    "listing": [
                        self.get_modified_output(x, listing_outdir)
                        for x in output["listing"]
                    ]
                }
            else:
                output = dict(
                    (key, self.get_modified_output(value, outdir))
                    for key, value in output
                )

        return output

    def create_results(self, outputs: dict, outdir: str) -> dict:
        return dict(
            (key, self.get_modified_output(value, outdir))
            for key, value in outputs.items()
            if not key.endswith("___dxfiles")
        )
