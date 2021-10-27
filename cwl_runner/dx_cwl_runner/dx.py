from datetime import datetime
import json
import os
from typing import Tuple

import dxpy

from dx_cwl_runner import utils
from dx_cwl_runner.utils import Log


def dx_file_to_path(f: dxpy.DXFile) -> str:
    desc = f.describe()
    return os.path.join(desc["folder"], desc["name"])


class Dx:
    def __init__(self, log: Log):
        self.log = log
        self._check_tools()
        self._cache = {}
        with open(
            f"{os.environ.get('HOME')}/.dnanexus_config/environment.json"
        ) as dx_env:
            env = json.load(dx_env)
            self.current_dx_folder = env.get("DX_CLI_WD", "/")
        self.current_dx_project = dxpy.PROJECT_CONTEXT_ID
        self.test_folder = os.path.join(
            self.current_dx_folder,
            f"cwl_runner_{datetime.now().strftime('%Y-%d-%m_%H-%M-%S')}",
        )
        self.inputs_folder = os.path.join(self.test_folder, "inputs")
        self.outputs_folder = os.path.join(self.test_folder, "outputs")
        if log.dryrun:
            log.info(f"creating inputs folder {self.inputs_folder}")
        else:
            proj = dxpy.DXProject(self.current_dx_project)
            proj.new_folder(self.inputs_folder, parents=True)

    def _check_tools(self):
        utils.check_tool("dx", "dx v([\\d.]+)", "dnanexus/dx-toolkit", self.log, lib="dxpy")

    def _get_cache(self, location: str) -> str:
        if location in self._cache:
            return self._cache[location]

    def _add_to_cache(self, local: str, remote: str):
        self._cache[local] = remote

    def get_file_name(self, file_id: str) -> str:
        dx_file = dxpy.DXFile(file_id, project=self.current_dx_project)
        return dx_file.describe(fields={"name": True}).get("name")

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
                folder=self.inputs_folder,
                properties={"checksum": checksum},
                describe=True,
            )
        )
        if len(existing) == 1:
            dx_file = existing[0]
        elif len(existing) > 1:
            raise Exception(
                f"found multiple files in {self.inputs_folder} with the same checksum: {existing}"
            )
        else:
            self.log.debug(f"uploading {file} to {self.inputs_folder}")
            if self.log.dryrun:
                dx_file = dxpy.DXFile(
                    "file-XXXXXXXXXXXXXXXXXXXXXXXX", self.current_dx_project
                )
            else:
                dx_file = dxpy.upload_local_file(
                    file, folder=self.inputs_folder, properties={"checksum": checksum}
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
        base_folder = os.path.join(self.inputs_folder, basename)
        dx_files = dict(
            (os.path.relpath(dx_file_to_path(f), base_folder), f)
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

    def download_file(self, file_id, output_path, project_id=None):
        dxpy.download_dxfile(
            file_id, output_path, project=project_id or self.current_dx_project
        )

    def download_folder(self, folder, outdir, project_id=None):
        dxpy.download_folder(
            project_id or self.current_dx_project, outdir, folder, overwrite=False
        )

    def run_executable(self, id, inputs: dict) -> Tuple[dict, str]:
        handler: dxpy.DXApplet = dxpy.get_handler(id)
        exe: dxpy.DXJob = handler.run(inputs, folder=self.outputs_folder)
        self.log.debug(f"Waiting for {exe.id} to finish...")
        exe.wait_on_done()
        desc = exe.describe()
        # TODO: is there an equivalent dxpy function?
        log = utils.run_cmd(f"dx watch {exe.id} --quiet", self.log)
        return desc, log
