import functools
import os
import tempfile
import hashlib
import time
import dxpy
from typing import List, Dict, Tuple, Set

from threading import Lock

from dxcint.Logger import Logger

logger = Logger.make(name=__name__, verbosity="info")

DEFAULT_INSTANCE_TYPE = "mem1_ssd1_v2_x4"


def rm_suffix(original_string: str, suffix: str) -> str:
    if not original_string:
        return ""
    if not suffix:
        return original_string
    if original_string[-1] == suffix[-1]:
        return rm_suffix(original_string[:-1], suffix[:-1])
    else:
        raise ValueError("rm_suffix(): suffix is not present in the original string")


def rm_prefix(original_string: str, prefix: str) -> str:
    if not original_string:
        return ""
    if not prefix:
        return original_string
    if original_string[0] == prefix[0]:
        return rm_prefix(original_string[1:], prefix[1:])
    else:
        raise ValueError("rm_suffix(): prefix is not present in the original string")


def async_retry(max_retries: int = 5, delay: int = 5):
    """
    A decorator function to perform async retry of a decorated method.
    Args:
        max_retries: int. Number of maximum retries.
        delay: int. Amount of time to sleep in seconds between retries.
    Returns: A result of a decorated callable.
    """

    def async_retry_inner(func):
        @functools.wraps(func)
        def async_retry_wrapper(*args, **kwargs):
            try:
                lock = args[
                    0
                ].context.lock  # For methods of a class with Context property
            except AttributeError:
                lock = Lock()  # For other functions.
            for i in range(0, max_retries):
                try:
                    with lock:
                        logger.info(
                            f"Retry {i} for function `{func.__name__}({args}, {kwargs})`"
                        )
                    ret_value = func(*args, **kwargs)
                    return ret_value
                except Exception as e:
                    with lock:
                        logger.info(
                            f"Error when running an async retry for function `{func.__name__}`\n"
                            f"Error CONTENT: {e}\n"
                            f"With ARGS: {args}\n"
                            f"With KWARGS: {kwargs}\n"
                            f"Retry in {delay} sec"
                        )
                    time.sleep(delay)
            else:
                raise Exception(
                    f"Failed after {max_retries} retries for function `{func.__name__}`\n"
                    f"With ARGS: {args}\n"
                    f"With KWARGS: {kwargs}"
                )

        return async_retry_wrapper

    return async_retry_inner


def list_dx_folder(project, folder, folder_cache):
    # get shallow listing of remote folder
    if isinstance(project, str):
        project = dxpy.DXProject(project)
    key = (project.get_id(), folder)
    if key in folder_cache:
        return folder_cache[key]
    contents = project.list_folder(folder)
    files: List[dict] = [
        {"$dnanexus_link": {"id": obj["id"], "project": project.get_id()}}
        for obj in contents["objects"]
        if obj["id"].startswith("file-")
    ]
    dirs: List[dict] = [
        {"type": "Folder", "uri": "dx://{}:{}".format(project.get_id(), folder)}
        for folder in contents["folders"]
    ]
    listing = files + dirs
    folder_cache[key] = listing
    return listing


def dict_compare(actual: Dict, expected: Dict) -> Tuple[Set, Set, Set, Set]:
    d1_keys = set(actual.keys())
    d2_keys = set(expected.keys())
    shared_keys = d1_keys.intersection(d2_keys)
    added = d1_keys - d2_keys
    removed = d2_keys - d1_keys
    modified = {}
    for o in shared_keys:
        if not isinstance(actual[o], type(expected[o])):
            modified[o] = (actual[o], expected[o])
            continue
        if isinstance(actual[o], dict):
            a, r, m, s = dict_compare(actual[o], expected[o])
            if not r and not m:  # expected dict cannot contain more keys, actual can
                continue
        elif actual[o] == expected[o]:
            continue
        modified[o] = (actual[o], expected[o])
    same = set(o for o in shared_keys if actual[o] == expected[o])
    return added, removed, modified, same


def sort_dicts(a: List[Dict], b: List[Dict]) -> Tuple[List[Dict], List[Dict]]:
    """Sort two lists of dicts to make them comparable. Given lists of dicts a and b:

    1. get the set of all keys in all dicts in both lists
    2. expand each dict into a list of tuples where the first element is the key and the
    second element is the value (which may be None if the dict does not contain the key)
    3. sort each list
    4. remove the tuples with None values from each list
    5. turn each list of tuples back into a list of dicts.
    """
    all_keys = list(
        sorted(
            set(k for x in a for k in x.keys()) | set(k for x in b for k in x.keys())
        )
    )
    return (
        [
            dict((k, v) for k, v in x if v is not None)
            for x in list(sorted([(k, i.get(k)) for k in all_keys] for i in a))
        ],
        [
            dict((k, v) for k, v in x if v is not None)
            for x in list(sorted([(k, j.get(k)) for k in all_keys] for j in b))
        ],
    )


def sort_maybe_mixed(seq):
    try:
        # may fail if the lists contain mutliple types of values
        return list(sorted(seq))
    except Exception:
        d = dict((str(x), x) for x in seq)
        sorted_keys = list(sorted(d.keys()))
        return [d[k] for k in sorted_keys]


def download_dxfile(dxfile, file_cache):
    key = (dxfile.get_proj_id(), dxfile.get_id())
    if key in file_cache:
        return file_cache[key]
    # the result is a file - download it and extract the contents
    dlpath = os.path.join(tempfile.mkdtemp(), dxfile.describe()["name"])
    dxpy.download_dxfile(dxfile, dlpath)
    try:
        with open(dlpath, "r") as inp:
            contents = str(inp.read()).strip()
            file_cache[key] = contents
            return contents
    finally:
        if os.path.exists(dlpath):
            os.remove(dlpath)


def link_to_dxfile(link, project_id):
    fields = link["$dnanexus_link"]
    if isinstance(fields, str):
        return dxpy.DXFile(fields, project_id)
    else:
        return dxpy.DXFile(fields["id"], fields.get("project", project_id))


def get_checksum(contents, algo):
    try:
        m = hashlib.new(algo)
        m.update(contents)
        checksum = m.digest()
        return f"{algo}${checksum}"
    except Exception:
        print("python does not support digest algorithm {}".format(algo))
        return None
