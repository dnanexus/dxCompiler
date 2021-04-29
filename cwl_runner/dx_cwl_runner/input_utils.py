import json
import os

from dx_cwl_runner.dx import Dx
from dx_cwl_runner.utils import Log


def upload_file(path: str, basedir: str, dx: Dx) -> str:
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


def upload_dir(path: str, basedir: str, dx: Dx) -> (str, str):
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


def get_modified_input(i, basedir: str, dx: Dx):
    if type(i) is list:
        return [get_modified_input(x, basedir, dx) for x in i]

    # File and Directory are always dicts. There may also be record types.
    # We differentiate based on the value of "class".
    if type(i) is dict:
        cls = i.get("class")
        if cls == "File":
            if "location" in i:
                i["location"] = upload_file(i["location"], basedir, dx)
            elif "path" in i:
                i["location"] = upload_file(i.pop("path"), basedir, dx)
            secondary_files = i.get("secondaryFiles")
            if secondary_files is not None:
                i["secondaryFiles"] = [get_modified_input(x, basedir, dx) for x in secondary_files]
        elif cls == "Directory":
            location = i.get("location", i.pop("path", None))
            if location is not None:
                location, basename = upload_dir(i["location"], basedir, dx)
                i["location"] = location
                if "basename" not in i and basename is not None:
                    i["basename"] = basename
            else:
                listing = i["listing"]
                if listing is None:
                    raise Exception(f"Directory is missing a 'location', 'path', or 'listing': {i}")
                i["listing"] = get_modified_input(listing, basedir, dx)
        else:
            dict((k, get_modified_input(v, basedir, dx)) for k, v in i.items())

    return i


def create_dx_input(jobfile: str, basedir: str, dx: Dx) -> (dict, str):
    # issues:
    #         if file does not exist, exception is thrown and no json is generated, even if some files were uploaded
    with open(jobfile) as input_file:
        inputs = json.load(input_file)
        process_file = inputs.pop("cwl:tool", None)
        return dict(
            (k, get_modified_input(v, basedir, dx))
            for k, v in inputs.items()
        ), process_file


def write_dx_input(new_input: dict, jobfile: str, log: Log) -> str:
    js = json.dumps(new_input, indent=4)
    if log.dryrun:
        log.log(f"writing input to {jobfile}:\n{js}")
    else:
        with open(jobfile, 'w') as dx_input:
            dx_input.write(js)
    return jobfile


def get_new_dx_input(input_path: str) -> str:
    basename = os.path.basename(input_path)
    new_input_filename = f"{os.path.splitext(basename)[0]}.dx.json"
    return new_input_filename
