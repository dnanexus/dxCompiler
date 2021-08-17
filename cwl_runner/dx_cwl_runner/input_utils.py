import json
import os
from posixpath import basename
import yaml
from typing import List, Optional, Tuple

from dx_cwl_runner.dx import Dx
from dx_cwl_runner.utils import Log, run_cmd


def upload_file(path: str, basedir: str, dx: Dx, cache: dict) -> str:
    if path in cache:
        return cache[path]
    if path.startswith("dx://"):
        return path
    if not os.path.isabs(path):
        path = os.path.join(basedir, path)
    # see if we've already found/uploaded this file
    if not os.path.exists(path):
        raise Exception(f"path does not exist: {path}")
    if not os.path.isfile(path):
        raise Exception(f"not a file {path}")
    return dx.find_or_upload_file(path)


def upload_dir(path: str, basedir: str, dx: Dx, cache: dict) -> Tuple[str, Optional[str]]:
    if path in cache:
        return cache[path]
    if path.startswith("dx://"):
        return path, None
    if not os.path.isabs(path):
        path = os.path.join(basedir, path)
    # see if we've already found/uploaded this file
    if not os.path.exists(path):
        raise Exception(f"path does not exist: {path}")
    if not os.path.isdir(path):
        raise Exception(f"not a file {path}")
    dx_uri = dx.find_or_upload_dir(path)
    basename = os.path.basename(path)
    return dx_uri, basename


def get_modified_input(i, basedir: str, dx: Dx, cache: dict):
    if type(i) is list:
        return [get_modified_input(x, basedir, dx, cache) for x in i]

    # File and Directory are always dicts. There may also be record types.
    # We differentiate based on the value of "class".
    if isinstance(i, dict):
        cls = i.get("class")
        if cls == "File":
            if "location" in i:
                i["location"] = upload_file(i["location"], basedir, dx, cache)
            elif "path" in i:
                i["location"] = upload_file(i.pop("path"), basedir, dx, cache)
            secondary_files = i.get("secondaryFiles")
            if secondary_files is not None:
                i["secondaryFiles"] = [
                    get_modified_input(x, basedir, dx, cache) for x in secondary_files
                ]
        elif cls == "Directory":
            location = i.get("location", i.pop("path", None))
            if location is not None:
                location, basename = upload_dir(i["location"], basedir, dx, cache)
                i["location"] = location
                if "basename" not in i and basename is not None:
                    i["basename"] = basename
            else:
                listing = i["listing"]
                if listing is None:
                    raise Exception(
                        f"Directory is missing a 'location', 'path', or 'listing': {i}"
                    )
                i["listing"] = get_modified_input(listing, basedir, dx, cache)
        else:
            i = dict((k, get_modified_input(v, basedir, dx, cache)) for k, v in i.items())

    return i


def create_dx_input(
    process_file: str, jobfile: str, basedir: str, tmpdir: str, dx: Dx, log: Log
) -> Tuple[str, str]:
    # issues:
    #   if file does not exist, exception is thrown and no json is generated, even if some
    #   files were uploaded
    
    # use a cache for uploaded files to make sure we upload each file/directory only once
    cache = {}

    # parse the inputs file
    with open(jobfile) as input_file:
        inputs = json.load(input_file)
        process_file = inputs.pop("cwl:tool", process_file)
        dx_inputs = dict(
            (k, get_modified_input(v, basedir, dx, cache))
            for k, v in inputs.items()
        )
        dx_input_file = os.path.join(tmpdir, os.path.basename(jobfile))
        write_dx_input(dx_inputs, dx_input_file, log)
    
    # check if the CWL has any hard-coded paths
    paths_referenced = get_paths_referenced(process_file, log)
    if paths_referenced:
        pathmap = dict(
            (
                p.get("location", p.get("path")), 
                get_modified_input(p, basedir, dx, cache).get("location")
            )
            for p in paths_referenced
            if "location" in p or "path" in p
        )

        if process_file.endswith(".json"):
            lib = json
        else:
            lib = yaml

        with open(process_file, "rt") as inp:
            cwl = lib.load(inp)
        dx_cwl = replace_paths(cwl, pathmap)
        process_file = os.path.join(tmpdir, os.path.basename(process_file))
        with open(process_file, "wt") as out:
            lib.write(dx_cwl, out)

    return (process_file, dx_input_file)


def write_dx_input(new_input: dict, jobfile: str, log: Log):
    js = json.dumps(new_input, indent=4)
    if log.dryrun:
        log.log(f"writing input to {jobfile}:\n{js}")
    else:
        with open(jobfile, "w") as dx_input:
            dx_input.write(js)


def get_new_dx_input(input_path: str) -> str:
    basename = os.path.basename(input_path)
    new_input_filename = f"{os.path.splitext(basename)[0]}.dx.json"
    return new_input_filename


def get_paths_referenced(process_file: str, log: Log) -> Optional[dict]:
    deps = json.loads(run_cmd(f"cwltool --print-deps {process_file}"))
    deps.get("secondaryFiles")


def replace_paths(obj, pathmap):
    if isinstance(obj, dict):
        if obj.get("class") in ("File", "Directory"):
            key = obj.get("location", obj.get("path"))
            if key in pathmap:
                obj["location"] = pathmap[key]
                return obj
        else:
            return dict((k, replace_paths(v)) for k, v in obj.items())
    elif isinstance(obj, list):
        return [replace_paths(item, pathmap) for item in obj]
    else:
        return obj
