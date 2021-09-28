from datetime import datetime
import glob
import json
import os
from pathlib import Path
try:
    from packaging.version import parse as parse_version
except:
    parse_version = None
import urllib.request

import dxpy

from dx_cwl_runner import utils


def get_latest_release():
    with urllib.request.urlopen(
            "https://api.github.com/repos/dnanexus/dxCompiler/releases/latest"
    ) as url:
        return json.loads(url.read().decode())


def download_dx_compiler(release_json) -> Path:
    asset = release_json["assets"][0]
    filename = Path(asset["name"])
    urllib.request.urlretrieve(asset["url"], filename)
    return filename


def pick_latest_version(jars) -> (Path, str):
    if parse_version:
        versions = dict(
            (parse_version(utils.run_cmd(f"java -jar {jar} version")), jar)
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


def get_compiler_jar(log) -> (Path, str):
    try:
        latest_release = get_latest_release()
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
            compiler_jar, compiler_version = pick_latest_version(jars)
            force_update = True
        elif latest_release:
            compiler_jar = download_dx_compiler(latest_release)
            compiler_version = latest_release["tag"]
            force_update = None
        else:
            raise Exception(
                "dxCompiler JAR not found in the 'DX_COMPILER_JAR' environment variable, nor "
                "in the current directory, and the latest version could not be retrieved from GitHub"
            )

    if force_update is not None and latest_release:
        latest_version_str = latest_release["tag_name"]
        if parse_version:
            current_version = parse_version(compiler_version)
            latest_version = parse_version(latest_version_str)
            if current_version < latest_version:
                if force_update:
                    compiler_jar = download_dx_compiler(latest_release)
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


class Dx:
    _compiler_jar = None

    @staticmethod
    def _get_compiler_jar(log) -> Path:
        if Dx._compiler_jar is None:
            Dx._compiler_jar = get_compiler_jar(log)
        return Dx._compiler_jar

    @staticmethod
    def dx_file_to_path(f: dxpy.DXFile) -> str:
        desc = f.describe()
        return os.path.join(desc["folder"], desc["name"])

    def __init__(self, config: utils.Log):
        self.log = config
        self.compiler_jar = Dx._get_compiler_jar(config)
        self.cache = {}
        with open(
            f"{os.environ.get('HOME')}/.dnanexus_config/environment.json"
        ) as dx_env:
            env = json.load(dx_env)
            self.current_dx_folder = env.get("DX_CLI_WD", "/")
        self.current_dx_project = dxpy.PROJECT_CONTEXT_ID
        self.test_folder = os.path.join(
            self.current_dx_folder,
            f"cwl_runner_{datetime.now().strftime('%Y.%d.%m_%H-%M-%S')}",
        )
        mk_folder_cmd = f"dx mkdir -p {self.current_dx_project}:{self.test_folder}"
        if config.dryrun:
            config.log(f"creating folder using command: {mk_folder_cmd}")
        else:
            utils.run_cmd(mk_folder_cmd, config.verbose)

    def _get_cache(self, location: str) -> str:
        if location in self.cache:
            return self.cache[location]

    def _add_to_cache(self, local: str, remote: str):
        self.cache[local] = remote

    def get_file_name(self, outdir: str, file_id: str) -> str:
        dx_file = dxpy.DXFile(file_id, project=self.current_dx_project)
        platform_file_name = dx_file.describe(fields={"name": True}).get("name")
        file_name = os.path.join(outdir, platform_file_name)
        counter = 0
        while os.path.exists(file_name):
            counter += 1
            file_name = f"{os.path.join(outdir, platform_file_name)}({counter})"
        return file_name

    def find_or_upload_file(self, file: str) -> str:
        # the file may have been resolved previously
        dx_uri = self._get_cache(file)
        if dx_uri is not None:
            return dx_uri

        # the file may have been uploaded in a previous run - search by md5
        checksum = utils.get_checksum(file)
        existing = list(
            dxpy.find_data_objects(
                classname="file",
                project=self.current_dx_project,
                folder=self.test_folder,
                properties={"checksum": checksum},
                describe=True,
            )
        )
        if len(existing) == 1:
            dx_file = existing[0]
        elif len(existing) > 1:
            raise Exception(
                f"found multiple files in {self.test_folder} with the same checksum: {existing}"
            )
        else:
            self.log.debug(f"uploading {file} to {self.test_folder}")
            if self.log.dryrun:
                dx_file = dxpy.DXFile(
                    "file-XXXXXXXXXXXXXXXXXXXXXXXX", self.current_dx_project
                )
            else:
                dx_file = dxpy.upload_local_file(
                    file, folder=self.test_folder, properties={"checksum": checksum}
                )
        dx_uri = f"dx://{dx_file.get_proj_id()}:{dx_file.get_id()}"
        self._add_to_cache(file, dx_uri)
        return dx_uri

    def find_or_upload_dir(self, dir: str) -> str:
        # the file may have been resolved previously
        dx_uri = self._get_cache(dir)
        if dx_uri is not None:
            return dx_uri

        local_files = dict(
            (os.path.relpath(p, dir), (p, utils.get_checksum(p)))
            for p in [
                os.path.join(dirpath, filename)
                for dirpath, dirnames, filenames in os.walk(dir)
                for filename in filenames
            ]
        )

        # the folder may have been uploaded in a previous run -
        # if so, make sure the contents are identical; if not,
        # create a new folder with a unique name
        basename = os.path.basename(dir)
        base_folder = os.path.join(self.test_folder, basename)
        dx_files = dict(
            (os.path.relpath(self.dx_file_to_path(f), base_folder), f)
            for f in dxpy.find_data_objects(
                classname="file",
                project=self.current_dx_project,
                folder=base_folder,
                describe={"name": True, "folder": True, "properties": True},
            )
        )

        if set(local_files.keys()) != set(dx_files.keys()):
            match = False
        else:
            # the paths match, but compare the checksums to make sure
            for relpath, (local_path, checksum) in local_files.items():
                dx_file = dx_files[relpath]
                dx_checksum = dx_file.describe().get("properties", {}).get("checksum")
                if dx_checksum is None or checksum != dx_checksum:
                    # there is not a perfect match, so upload all files to a new folder
                    match = False
                    break
            else:
                # all paths match
                match = True

        if not match:
            # create a new folder
            proj = dxpy.DXProject(self.current_dx_project)
            i = 0
            while True:
                base_folder = f"{basename}_{i}"
                # check if the folder exists by listing it and seeing if it contains
                # any files or folders
                list_obj = proj.list_folder(base_folder)
                if not (
                    list_obj and (list_obj.get("objects") or list_obj.get("folders"))
                ):
                    i += 1
                    break

            self.log.debug(f"uploading contents of {dir} to {base_folder}")

            for relpath, (local_path, checksum) in local_files.items():
                folder = os.path.join(base_folder, os.path.dirname(relpath))
                self.log.debug(f"uploading {local_path} to {folder}")
                if not self.log.dryrun:
                    dxpy.upload_local_file(
                        local_path, folder=folder, properties={"checksum": checksum}
                    )

        return f"dx://{self.current_dx_project}:{base_folder}"

    def check_outdir(self, dir: str, create: bool = True):
        if not os.path.exists(dir):
            self.log.log(f"Creating output directory {dir}")
            if create:
                os.mkdir(dir)
        elif not os.path.isdir(dir):
            raise Exception(f"Outdir {dir} is not a directory!")
        elif not os.access(dir, os.W_OK):
            raise Exception(f"You need write access to the outdir repository ({dir})!")
