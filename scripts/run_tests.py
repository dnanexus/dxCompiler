#!/usr/bin/env python3
import argparse
from collections import namedtuple
import dxpy
import glob
import hashlib
import json
import os
import re
import sys
import subprocess
import tempfile
from termcolor import cprint
import time
import traceback
import yaml
from dxpy.exceptions import DXJobFailureError

import util

here = os.path.dirname(sys.argv[0])
top_dir = os.path.dirname(os.path.abspath(here))
test_dir = os.path.join(os.path.abspath(top_dir), "test")

git_revision = subprocess.check_output(
    ["git", "describe", "--always", "--dirty", "--tags"]
).strip()
test_files = {}
test_failing = {
    "bad_status",
    "bad_status2",
    "just_fail_wf",
    "missing_output",
    "docker_retry",
    "argument_list_too_long",
}

wdl_v1_list = [
    # calling native dx applets/apps
    "call_native_v1",
    "call_native_app",
    "cast",
    "dict",
    "instance_types",
    "linear_no_expressions",
    "linear",
    "optionals",
    "optionals3",
    "spaces_in_file_paths",
    "strings",
    "runtime_vs_static_type",
    "wf_person",
    "call_level2",
    "environment_passing_deep_nesting",
    "optional_output",
    "unpassed_default_arg",
    # workflows with nested blocks
    "two_levels",
    "three_levels",
    "four_levels",
    "param_passing",
    "nested_scatter",
    # Array input with no values
    "empty_array",
    # Map with a File key
    "map_file_key",
    # defaults and parameter passing
    "top",
    "subworkflow_with_default",
    # can we download from a container?
    "download_from_container",
    # input file with pairs
    "echo_pairs",
    "array_structs",
    # Missing optional output files, returned as none, instead
    # of an error
    "missing_optional_output_file",
    # calling with an optional argument not specified
    "scatter_subworkflow_with_optional",
    # streaming
    "streaming_inputs",
    # input/output linear_no_expressions
    "wf_with_input_expressions",
    "wf_with_output_expressions",
    # bug regression tests
    "nested_pairs",  # APPS-370
    "apps_378",
    "apps_384",
    "diff_stream_and_download",  # APPS-288
    "apps_573",
    "apps_612",
    "nested_optional",
    "struct_deref",  # APPS-615

    # manifests
    "simple_manifest",
    "complex_manifest",
    "view_and_count_manifest",

    # workflow with output files created by expressions
    "upload_workflow_files"
]

wdl_v1_1_list = [
    "v1_1_dict",

    # bug regression tests
    "apps_579_boolean_flag_expr",
    "apps_579_string_substitution_expr"
]

# docker image tests
docker_test_list = [
    "broad_genomics",
    "biocontainers",
    "private_registry",
    "native_docker_file_image",
    "native_docker_file_image_gzip",
    "samtools_count",
    "dynamic_docker_image",
    "ecr_docker",
]

# wdl draft-2
draft2_test_list = [
    "advanced",
    "bad_status",
    "bad_status2",
    "just_fail_wf",
    "call_with_defaults1",
    "call_with_defaults2",
    "conditionals_base",
    "files",
    "files_with_the_same_name",
    "hello",
    "shapes",
    # this test cannot be enabled yet, because we
    # don't yet support overriding task inputs
    #"population",

    # multiple library imports in one WDL workflow
    "multiple_imports",
    # subworkflows
    "conditionals2",
    "modulo",
    "movies",
    "subblocks2",
    "subblocks",
    "var_type_change",
    "outer",
    # calling native dx applets/apps
    # We currently do not have a code generator for draft-2, so cannot import dx_extern.wdl.
    # "call_native",
    "write_lines_bug",
]

single_tasks_list = [
    "add3",
    "diff2files",
    "empty_stdout",
    "sort_file",
    "symlinks_wc",
    "DiskSpace2",
    "echo_line_split",
    "opt_array",
    "stream_diff_v1",
]

cwl_tools = [
    "cat",  # hello world tool
    "tar_files",
]

cwl_conformance_tools = [
    os.path.basename(path)[:-4]
    for path in glob.glob(os.path.join(test_dir, "cwl_conformance", "tools", "*.cwl"))
]
cwl_conformance_workflows = [
    os.path.basename(path)[:-9]
    for path in glob.glob(os.path.join(test_dir, "cwl_conformance", "workflows_packed", "*.cwl.json"))
]
cwl_directory_tools = [
    "dir",
    "dir2",
    "dir3",
    "dir4",
    "dir5",
    "dir6"
]

# Tests run in continuous integration. We remove the native app test,
# because we don't want to give permissions for creating platform apps.
ci_test_list = [
    # WDL tests
    "advanced",
    # We currently do not have a code generator for draft-2, so cannot import dx_extern.wdl.
    # "call_native",
    "call_with_defaults1",
    "trains",
    "files",
    # CWL tests
    "cat",
]

special_flags_list = [
    "add2",  # test the ignoreReuse flag
    "add_many",  # tests the delayWorkspaceDestruction flag
    "inc_range",  # check that runtime call to job/analysis pass the delayWorkspaceDestruction flag
]

# these are the examples from the documentation
doc_tests_list = [
    "bwa_mem"
]

medium_test_list = (
    wdl_v1_list + wdl_v1_1_list + docker_test_list + special_flags_list + cwl_tools
)
large_test_list = (
    medium_test_list + draft2_test_list + single_tasks_list + doc_tests_list
)

test_suites = {
    "CI": ci_test_list,
    "M": medium_test_list,
    "L": large_test_list,
    "tasks": single_tasks_list,
    "draft2": draft2_test_list,
    "docker": docker_test_list,
    "native": ["call_native", "call_native_v1"],
    "docs": doc_tests_list,
    "cwl_conformance": cwl_conformance_tools + cwl_conformance_workflows,
    "cwl_dir": cwl_directory_tools
}

# Tests with the reorg flags
test_reorg = {
    "dict",
    "strings",
    "test_reorg",
    "test_reorg_no_config"
}
test_defaults = set()
test_unlocked = {
    "array_structs",
    "cast",
    "call_with_defaults1",
    "files",
    "hello",
    "path_not_taken",
    "optionals",
    "shapes",
    "population"
}
test_project_wide_reuse = {
    "add2",
    "add_many"
}
test_separate_outputs = {
    "localization"
}

test_import_dirs = ["A"]
TestMetaData = namedtuple("TestMetaData", ["name", "kind"])
TestDesc = namedtuple(
    "TestDesc",
    ["name", "kind", "source_file", "raw_input", "dx_input", "results", "extras"],
)

# Test with -waitOnUpload flag
test_upload_wait = {
    "upload_wait"
}


def read_json_file(path):
    with open(path, "r") as fd:
        data = fd.read()
        d = json.loads(data)
        return d


def verify_json_file(path):
    try:
        read_json_file(path)
    except:
        raise RuntimeError("Error verifying JSON file {}".format(path))


# Search a WDL file with a python regular expression.
# Note this is not 100% accurate.
#
# Look for all tasks and workflows. If there is exactly
# one workflow, this is a WORKFLOW. If there are no
# workflows and exactly one task, this is an APPLET.
task_pattern_re = re.compile(r"^(task)(\s+)(\w+)(\s+){")
wf_pattern_re = re.compile(r"^(workflow)(\s+)(\w+)(\s+){")


def get_wdl_metadata(filename):
    workflows = []
    tasks = []
    with open(filename, "r") as fd:
        for line in fd:
            m = re.match(wf_pattern_re, line)
            if m is not None:
                workflows.append(m.group(3))
            m = re.match(task_pattern_re, line)
            if m is not None:
                tasks.append(m.group(3))
    if len(workflows) > 1:
        raise RuntimeError("WDL file {} has multiple workflows".format(filename))
    if len(workflows) == 1:
        return TestMetaData(name=workflows[0], kind="workflow")
    assert len(workflows) == 0
    if len(tasks) == 1:
        return TestMetaData(name=tasks[0], kind="applet")
    if os.path.basename(filename).startswith("library") or os.path.basename(
        filename
    ).endswith("_extern"):
        return
    raise RuntimeError(
        "{} is not a valid WDL test, #tasks={}".format(filename, len(tasks))
    )


def get_cwl_metadata(filename, tname):
    with open(filename, "r") as fd:
        doc = yaml.safe_load(fd)

    if doc["class"] == "CommandLineTool":
        name = doc.get("id", tname)
        return TestMetaData(name=name, kind="applet")

    raise RuntimeError("{} is not a valid CWL tool test".format(filename))


def get_cwl_json_metadata(filename, tname):
    with open(filename, 'r') as fd:
        doc = json.load(fd)

    if "class" in doc:
        if doc["class"] == "Workflow":
            # the workflow id in a packed CWL file is always "main"
            # so we use the test name instead
            return TestMetaData(name=tname, kind="workflow")
    elif "$graph" in doc:
        for proc in doc["$graph"]:
            if proc["id"] == "#main":
                if proc["class"] == "Workflow":
                    return TestMetaData(name=tname, kind="workflow")
                else:
                    break

    raise RuntimeError("{} is not a valid CWL workflow test".format(filename))


# Register a test name, find its inputs and expected results files.
def register_test(dir_path, tname, ext):
    global test_files
    if tname in test_suites.keys():
        raise RuntimeError(
            "Test name {} is already used by a test-suite, it is reserved".format(tname)
        )
    source_file = os.path.join(dir_path, tname + ext)
    if not os.path.exists(source_file):
        raise RuntimeError("Test file {} does not exist".format(source_file))
    if ext == ".wdl":
        metadata = get_wdl_metadata(source_file)
    elif ext == ".cwl":
        metadata = get_cwl_metadata(source_file, tname)
    elif ext == ".cwl.json":
        metadata = get_cwl_json_metadata(source_file, tname)
    else:
        raise RuntimeError("unsupported file type {}".format(ext))
    desc = TestDesc(
        name=metadata.name,
        kind=metadata.kind,
        source_file=source_file,
        raw_input=[],
        dx_input=[],
        results=[],
        extras=None,
    )

    # Verify the input file, and add it (if it exists)
    test_input = os.path.join(dir_path, tname + "_input.json")
    if os.path.exists(test_input):
        verify_json_file(test_input)
        desc.raw_input.append(test_input)
        desc.dx_input.append(os.path.join(dir_path, tname + "_input.dx.json"))
        desc.results.append(os.path.join(dir_path, tname + "_results.json"))

    # check if the alternate naming scheme is used for tests with multiple inputs
    i = 1
    while True:
        test_input = os.path.join(dir_path, tname + "_input{}.json".format(i))
        if os.path.exists(test_input):
            verify_json_file(test_input)
            desc.raw_input.append(test_input)
            desc.dx_input.append(
                os.path.join(dir_path, tname + "_input{}.dx.json".format(i))
            )
            desc.results.append(
                os.path.join(dir_path, tname + "_results{}.json".format(i))
            )
            i += 1
        else:
            break

    # Add an extras file (if it exists)
    extras = os.path.join(dir_path, tname + "_extras.json")
    if os.path.exists(extras):
        desc = desc._replace(extras=extras)

    test_files[tname] = desc
    return desc


######################################################################

# Same as above, however, if a file is empty, return an empty dictionary
def read_json_file_maybe_empty(path):
    if not os.path.exists(path):
        return {}
    else:
        return read_json_file(path)


def find_test_from_exec(exec_obj):
    dx_desc = exec_obj.describe()
    exec_name = dx_desc["name"].split(" ")[0]
    for tname, desc in test_files.items():
        if desc.name == exec_name:
            return tname
    raise RuntimeError("Test for {} {} not found".format(exec_obj, exec_name))


def link_to_dxfile(link, project):
    fields = link["$dnanexus_link"]
    if isinstance(fields, str):
        return dxpy.DXFile(fields, project)
    else:
        return dxpy.DXFile(fields["id"], fields.get("project", project))


file_cache = {}


def download_dxfile(dxfile):
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


# result may be a string with contents of the file, a dx link object, a
#   serialized CwlFile (an object with class="File") or a serialized
#   VFile (an object with type="File")
# expected_val is a serialized CwlFile
def compare_result_file(result, expected_val, field_name, tname, project, verbose=True):
    if "checksum" in expected_val:
        algo, expected_checksum = expected_val["checksum"].split("$")
    else:
        algo = None
        expected_checksum = None

    location = None
    size = None
    checksum = None
    secondary_files = None

    if isinstance(result, str):
        contents = result
    elif result.get("class", result.get("type")) == "File":
        contents = result.get("contents")
        location = result.get("location", result.get("path", result.get("uri")))
        size = result.get("size")
        checksum = result.get("checksum")
        secondary_files = set(result.get("secondaryFiles", []))
    elif "$dnanexus_link" in result:
        dxfile = link_to_dxfile(result, project)
        contents = download_dxfile(dxfile)
        location = os.path.join(dxfile.describe()["folder"], dxfile.describe()["name"])
    else:
        raise Exception("unsupported file value {}".format(result))

    expected_location = expected_val.get("location", expected_val.get("path"))
    expected_basename = expected_val.get("basename", os.path.basename(expected_location) if expected_location else None)
    if expected_basename:
        basename = os.path.basename(location) if location else None
        if basename != expected_basename:
            if verbose:
                cprint("Analysis {} gave unexpected results".format(tname), "red")
                cprint(
                    "Field {} should have location/path with basename ({}), actual = ({})".format(
                        field_name, expected_basename, basename
                    ),
                    "red"
                )
            return False

    # the result is a cwl File - match the contents, checksum, and/or size
    if "contents" in expected_val and contents != expected_val["contents"]:
        if verbose:
            cprint("Analysis {} gave unexpected results".format(tname), "red")
            cprint(
                "Field {} should have contents ({}), actual = ({})".format(
                    field_name, expected_val["contents"], result.get("contents")
                ),
                "red"
            )
        return False

    if "size" in expected_val:
        size |= (len(contents) if contents else None)
        if size != expected_val["size"]:
            if verbose:
                cprint("Analysis {} gave unexpected results".format(tname), "red")
                cprint(
                    "Field {} should have size ({}), actual = ({})".format(
                        field_name, expected_val["size"], len(contents)
                    ),
                    "red"
                )
            return False

    if expected_checksum:
        checksum |= (get_checksum(contents, algo) if contents else None)
        if checksum != expected_checksum:
            if verbose:
                cprint("Analysis {} gave unexpected results".format(tname), "red")
                cprint(
                    "Field {} should have checksum ({}), actual = ({})".format(
                        field_name, expected_checksum, checksum
                    ),
                    "red"
                )
            return False

    expected_secondary_files = expected_val.get("secondaryFiles")
    if expected_secondary_files:
        seconary_files_len = len(secondary_files) if secondary_files else 0
        if len(expected_secondary_files) == seconary_files_len:
            if verbose:
                cprint("Analysis {} gave unexpected results".format(tname), "red")
                cprint(
                    "Field {} should have secondaryFiles ({}), actual = ({})".format(
                        field_name, expected_secondary_files, secondary_files
                    ),
                    "red"
                )
            return False
        for expected in expected_secondary_files:
            for actual in secondary_files:
                if compare_result_path(
                        actual, expected, "{}.secondaryFiles".format(field_name), tname, project, verbose=False
                ):
                    secondary_files.remove(actual)
                    break
            else:
                if verbose:
                    cprint("Analysis {} gave unexpected results".format(tname), "red")
                    cprint(
                        "Field {} is missing secondaryFile ({}) from ({})".format(
                            field_name, expected, secondary_files
                        ),
                        "red"
                    )
                return False

    return True


folder_cache = {}


def list_dx_folder(project, folder):
    # get shallow listing of remote folder
    if isinstance(project, str):
        project = dxpy.DXProject(project)
    key = (project.get_id(), folder)
    if key in folder_cache:
        return folder_cache[key]
    contents = project.list_folder(folder)
    files = [
        {"$dnanexus_link": {"id": obj["id"], "project": project.get_id()}}
        for obj in contents["objects"]
        if obj["id"].startswith("file-")
    ]
    dirs = [
        {"type": "Folder", "uri": "dx://{}:{}".format(project.get_id(), folder)}
        for folder in contents["folders"]
    ]
    listing = files + dirs
    folder_cache[key] = listing
    return listing


# expected_val is a serialized CwlDirectory (an object with class="Directory")
# result may be a serialized CwlDirectory or a serialized VFolder (an object
#   with type="Folder")
def compare_result_directory(result, expected_val, field_name, tname, project, verbose=True):
    location = result.get("location", result.get("path", result.get("uri")))
    basename = None
    if location is not None and location.startswith("dx://"):
        project, folder = location[5:].split(":")
    else:
        folder = location
    if folder:
        basename = os.path.basename(folder)

    expected_location = expected_val.get("location", expected_val.get("path"))
    expected_basename = expected_val.get("basename", os.path.basename(expected_location) if expected_location else None)
    if expected_basename:
        if basename != expected_basename:
            if verbose:
                cprint("Analysis {} gave unexpected results".format(tname), "red")
                cprint(
                    "Field {} should have location/path with basename ({}), actual = ({})".format(
                        field_name, expected_basename, basename
                    ),
                    "red"
                )
            return False

    expected_listing = expected_val.get("listing")
    if expected_listing:
        if "listing" in result:
            listing = result["listing"]
        elif not folder:
            if verbose:
                cprint("Analysis {} gave unexpected results".format(tname), "red")
                cprint("Field {} is missing a folder, actual = ({})".format(field_name, result), "red")
            return False
        else:
            listing = list_dx_folder(project, folder)
        listing_len = len(listing) if listing else 0
        if len(expected_listing) != listing_len:
            if verbose:
                cprint("Analysis {} gave unexpected results".format(tname), "red")
                cprint(
                    "Field {} should have listing ({}), actual = ({})".format(
                        field_name, expected_listing, listing
                    ),
                    "red"
                )
            return False
        for expected in expected_listing:
            for actual in listing:
                if compare_result_path(
                        actual, expected, "{}.listing".format(field_name), tname, project, verbose=False
                ):
                    listing.remove(actual)
                    break
            else:
                if verbose:
                    cprint("Analysis {} gave unexpected results".format(tname), "red")
                    cprint(
                        "Field {} is missing item ({}) from listing ({})".format(
                            field_name, expected, listing
                        ),
                        "red"
                    )
                return False

    return True


def compare_result_path(result, expected_val, field_name, tname, project, verbose=True):
    cls = result.get("class", result.get("type"))
    expected_cls = result.get("class", result.get("type"))
    if cls == "File" and expected_cls == "File":
        return compare_result_file(result, expected_val, field_name, tname, project, verbose)
    elif cls in {"Directory", "Folder"} and expected_cls in {"Directory", "Folder"}:
        return compare_result_directory(result, expected_val, field_name, tname, project, verbose)
    else:
        return False


# Check that a workflow returned the expected result for
# a [key]
def validate_result(tname, exec_outputs: dict, key, expected_val, project):
    desc = test_files[tname]
    # Extract the key. For example, for workflow "math" returning
    # output "count":
    #    'math.count' -> count
    exec_name = key.split(".")[0]
    field_name_parts = key.split(".")[1:]

    field_name1 = ".".join(field_name_parts)
    # convert dots to ___
    field_name2 = "___".join(field_name_parts)
    if exec_name != tname:
        raise RuntimeError(
            "Key {} is invalid, must start with {} name".format(key, desc.kind)
        )
    try:
        # get the actual results
        if field_name1 in exec_outputs:
            result = exec_outputs[field_name1]
        elif field_name2 in exec_outputs:
            result = exec_outputs[field_name2]
        else:
            cprint(
                "field {} missing from executable results {}".format(
                    field_name1, exec_outputs
                ),
                "red",
            )
            return False
        if isinstance(result, dict) and "___" in result:
            result = result["___"]
        if isinstance(result, list) and isinstance(expected_val, list):
            result.sort()
            expected_val.sort()
        if isinstance(expected_val, dict) and (
            expected_val.get("class") in {"File", "Directory"} or
            expected_val.get("type") in {"File", "Folder"}
        ):
            compare_result_path(result, expected_val, field_name1, tname, project)
        else:
            if isinstance(result, dict) and result.get("type") == "File" and "uri" in result:
                result = result["uri"]
            if isinstance(result, dict) and "$dnanexus_link" in result:
                result = download_dxfile(link_to_dxfile(result, project.id))
            if str(result).strip() != str(expected_val).strip():
                cprint("Analysis {} gave unexpected results".format(tname), "red")
                cprint(
                    "Field {} should be ({}), actual = ({})".format(
                        field_name1, expected_val, result
                    ),
                    "red",
                )
                return False
        return True
    except Exception as e:
        print("exception message={}".format(e))
        return False


def get_checksum(contents, algo):
    try:
        m = hashlib.new(algo)
        m.update(contents)
        return m.digest()
    except:
        print("python does not support digest algorithm {}".format(algo))
        return None


def lookup_dataobj(tname, project, folder):
    desc = test_files[tname]
    wfgen = dxpy.bindings.search.find_data_objects(
        classname=desc.kind,
        name=desc.name,
        folder=folder,
        project=project.get_id(),
        limit=1,
    )
    objs = [item for item in wfgen]
    if len(objs) > 0:
        return objs[0]["id"]
    return None


# Build a workflow.
#
# wf             workflow name
# classpath      java classpath needed for running compilation
# folder         destination folder on the platform
def build_test(tname, project, folder, version_id, compiler_flags):
    desc = test_files[tname]
    print("build {} {}".format(desc.kind, desc.name))
    print("Compiling {} to a {}".format(desc.source_file, desc.kind))
    cmdline = [
        "java",
        "-jar",
        os.path.join(top_dir, "dxCompiler-{}.jar".format(version_id)),
        "compile",
        desc.source_file,
        "-force",
        "-folder",
        folder,
        "-project",
        project.get_id(),
    ]
    if "manifest" in desc.source_file:
        cmdline.append("-useManifests")
    cmdline += compiler_flags
    print(" ".join(cmdline))
    oid = subprocess.check_output(cmdline).strip()
    return oid.decode("ascii")


def ensure_dir(path):
    print("making sure that {} exists".format(path))
    if not os.path.exists(path):
        os.makedirs(path)


def wait_for_completion(test_exec_objs):
    print("awaiting completion ...")
    failures = []
    for exec_obj in test_exec_objs:
        tname = find_test_from_exec(exec_obj)
        desc = test_files[tname]
        try:
            exec_obj.wait_on_done()
            print("Executable {} succeeded".format(desc.name))
        except DXJobFailureError:
            if tname in test_failing:
                print("Executable {} failed as expected".format(desc.name))
            else:
                cprint("Error: executable {} failed".format(desc.name), "red")
                failures.append(tname)
    print("tools execution completed")
    return failures


# Run [workflow] on several inputs, return the analysis ID.
def run_executable(
    project, test_folder, tname, oid, debug_flag, delay_workspace_destruction
):
    desc = test_files[tname]

    def once(i):
        try:
            if tname in test_defaults or i < 0:
                print("  with empty input")
                inputs = {}
            else:
                print("  with input file: {}".format(desc.dx_input[i]))
                inputs = read_json_file(desc.dx_input[i])
            project.new_folder(test_folder, parents=True)
            if desc.kind == "workflow":
                exec_obj = dxpy.DXWorkflow(project=project.get_id(), dxid=oid)
            elif desc.kind == "applet":
                exec_obj = dxpy.DXApplet(project=project.get_id(), dxid=oid)
            else:
                raise RuntimeError("Unknown kind {}".format(desc.kind))

            run_kwargs = {}
            if debug_flag:
                run_kwargs = {
                    "debug": {
                        "debugOn": ["AppError", "AppInternalError", "ExecutionError"]
                    },
                    "allow_ssh": ["*"],
                }
            if delay_workspace_destruction:
                run_kwargs["delay_workspace_destruction"] = True

            return exec_obj.run(
                inputs,
                project=project.get_id(),
                folder=test_folder,
                name="{} {}".format(desc.name, git_revision),
                instance_type="mem1_ssd1_x4",
                **run_kwargs,
            )
        except Exception as e:
            print("exception message={}".format(e))
            return None

    def run(i):
        for _ in range(1, 5):
            retval = once(i)
            if retval is not None:
                return retval
            print("Sleeping for 5 seconds before trying again")
            time.sleep(5)
        else:
            raise RuntimeError("running workflow")

    n = len(desc.dx_input)
    if n == 0:
        return [run(-1)]
    else:
        return [run(i) for i in range(n)]


def extract_outputs(tname, exec_obj) -> dict:
    desc = test_files[tname]
    if desc.kind == "workflow":
        locked = tname not in test_unlocked
        if locked:
            return exec_obj["output"]
        else:
            stages = exec_obj["stages"]
            for snum in range(len(stages)):
                crnt = stages[snum]
                if crnt["id"] == "stage-outputs":
                    return stages[snum]["execution"]["output"]
            raise RuntimeError(
                "Analysis for test {} does not have stage 'outputs'".format(tname)
            )
    elif desc.kind == "applet":
        return exec_obj["output"]
    else:
        raise RuntimeError("Unknown kind {}".format(desc.kind))


def run_test_subset(
    project, runnable, test_folder, debug_flag, delay_workspace_destruction
):
    # Run the workflows
    test_exec_objs = []
    for tname, oid in runnable.items():
        desc = test_files[tname]
        print("Running {} {} {}".format(desc.kind, desc.name, oid))
        anl = run_executable(
            project, test_folder, tname, oid, debug_flag, delay_workspace_destruction
        )
        test_exec_objs.extend(anl)
    print("executables: " + ", ".join([a.get_id() for a in test_exec_objs]))

    # Wait for completion
    failed_execution = wait_for_completion(test_exec_objs)

    print("Verifying results")

    def verify_test(exec_obj, i):
        exec_desc = exec_obj.describe()
        tname = find_test_from_exec(exec_obj)
        if tname in test_failing:
            return None
        test_desc = test_files[tname]
        exec_outputs = extract_outputs(tname, exec_desc)
        if len(test_desc.results) > i:
            shouldbe = read_json_file_maybe_empty(test_desc.results[i])
            correct = True
            print("Checking results for workflow {} job {}".format(test_desc.name, i))
            for key, expected_val in shouldbe.items():
                correct = validate_result(tname, exec_outputs, key, expected_val, project)
            anl_name = "{}.{}".format(tname, i)
            if correct:
                print("Analysis {} passed".format(anl_name))
                return None
            else:
                return anl_name

    failed_verification = []
    for i, exec_obj in enumerate(test_exec_objs):
        failed_name = verify_test(exec_obj, i)
        if failed_name is not None:
            failed_verification.append(failed_name)

    if failed_execution or failed_verification:
        print("-----------------------------")
        if failed_execution:
            fexec = "\n".join(failed_execution)
            print(f"Tools failed execution: {len(failed_execution)}:\n{fexec}")
        if failed_verification:
            fveri = "\n".join(failed_verification)
            print(
                f"Tools failed results verification: {len(failed_verification)}:\n{fveri}"
            )
        # write failed tests to a file so we can easily re-run them next time
        with open(".failed", "wt") as out:
            out.write("\n".join(failed_execution + failed_verification))
        raise RuntimeError("Failed")


def print_test_list():
    test_list = "\n  ".join(sorted(key for key in test_files.keys()))
    print("List of tests:\n  {}".format(test_list))


# Choose set set of tests to run
def choose_tests(name):
    if name in test_suites.keys():
        return test_suites[name]
    if name == "All":
        return test_files.keys()
    if name in test_files.keys():
        return [name]
    # Last chance: check if the name is a prefix.
    # Accept it if there is exactly a single match.
    matches = [key for key in test_files.keys() if key.startswith(name)]
    if len(matches) > 1:
        raise RuntimeError(
            "Too many matches for test prefix {} -> {}".format(name, matches)
        )
    if len(matches) == 0:
        raise RuntimeError("Test prefix {} is unknown".format(name))
    return matches


# Find all the WDL test files, these are located in the 'test'
# directory. A test file must have some support files.
def register_all_tests(verbose: bool) -> None:
    for root, dirs, files in os.walk(test_dir):
        if os.path.basename(root).endswith("_ignore") or os.path.basename(root).endswith("_notimplemented"):
            continue
        for t_file in files:
            if t_file.endswith(".wdl") or t_file.endswith(".cwl"):
                base = os.path.basename(t_file)
                (fname, ext) = os.path.splitext(base)
            elif t_file.endswith(".cwl.json"):
                base = os.path.basename(t_file)
                fname = base[:-9]
                ext = ".cwl.json"
            else:
                continue

            if fname.startswith("library_"):
                continue
            if fname.endswith("_extern"):
                continue
            try:
                register_test(root, fname, ext)
            except Exception as e:
                if verbose:
                    print("Skipping file {} error={}".format(fname, e))


# Some compiler flags are test specific
def compiler_per_test_flags(tname):
    flags = []
    desc = test_files[tname]
    if tname not in test_unlocked:
        flags.append("-locked")
    if tname in test_reorg:
        flags.append("-reorg")
    if tname in test_project_wide_reuse:
        flags.append("-projectWideReuse")
    if tname in test_separate_outputs:
        flags.append("-separateOutputs")
    if tname in test_defaults and len(desc.raw_input) > 0:
        flags.append("-defaults")
        flags.append(desc.raw_input[0])
    if tname in test_upload_wait:
        flags.append("-waitOnUpload")
    else:
        for i in desc.raw_input:
            flags.append("-inputs")
            flags.append(i)
    if desc.extras is not None:
        flags += ["--extras", os.path.join(top_dir, desc.extras)]
    if tname in test_import_dirs:
        flags += ["--imports", os.path.join(top_dir, "test/imports/lib")]
    return flags


# Which project to use for a test
# def project_for_test(tname):

######################################################################


def native_call_dxni(project, applet_folder, version_id, verbose: bool):
    # build WDL wrapper tasks in test/dx_extern.wdl
    cmdline_common = [
        "java",
        "-jar",
        os.path.join(top_dir, "dxCompiler-{}.jar".format(version_id)),
        "dxni",
        "-force",
        "-folder",
        applet_folder,
        "-project",
        project.get_id(),
    ]
    if verbose:
        cmdline_common.append("--verbose")

    # draft-2 is not currently supported
    #     cmdline_draft2 = cmdline_common + [ "--language", "wdl_draft2",
    #                                         "--output", os.path.join(top_dir, "test/draft2/dx_extern.wdl")]
    #     print(" ".join(cmdline_draft2))
    #     subprocess.check_output(cmdline_draft2)

    cmdline_v1 = cmdline_common + [
        "-language",
        "wdl_v1.0",
        "-output",
        os.path.join(top_dir, "test/wdl_1_0/dx_extern.wdl"),
    ]
    print(" ".join(cmdline_v1))
    subprocess.check_output(cmdline_v1)


def dxni_call_with_path(project, path, version_id, verbose):
    # build WDL wrapper tasks in test/dx_extern.wdl
    cmdline = [
        "java",
        "-jar",
        os.path.join(top_dir, "dxCompiler-{}.jar".format(version_id)),
        "dxni",
        "-force",
        "-path",
        path,
        "-language",
        "wdl_v1.0",
        "-output",
        os.path.join(top_dir, "test/wdl_1_0/dx_extern_one.wdl"),
    ]
    if project is not None:
        cmdline.extend(["-project", project.get_id()])
    if verbose:
        cmdline.append("-verbose")
    print(" ".join(cmdline))
    subprocess.check_output(cmdline)


# Set up the native calling tests
def native_call_setup(project, applet_folder, version_id, verbose):
    native_applets = [
        "native_concat",
        "native_diff",
        "native_mk_list",
        "native_sum",
        "native_sum_012",
    ]

    # build the native applets, only if they do not exist
    for napl in native_applets:
        applet = list(
            dxpy.bindings.search.find_data_objects(
                classname="applet",
                name=napl,
                folder=applet_folder,
                project=project.get_id(),
            )
        )
        if len(applet) == 0:
            cmdline = [
                "dx",
                "build",
                os.path.join(top_dir, "test/applets/{}".format(napl)),
                "--destination",
                (project.get_id() + ":" + applet_folder + "/"),
            ]
            print(" ".join(cmdline))
            subprocess.check_output(cmdline)

    dxni_call_with_path(project, applet_folder + "/native_concat", version_id, verbose)
    native_call_dxni(project, applet_folder, version_id, verbose)

    # check if providing an applet-id in the path argument works
    first_applet = native_applets[0]
    results = dxpy.bindings.search.find_one_data_object(
        classname="applet",
        name=first_applet,
        folder=applet_folder,
        project=project.get_id(),
    )
    if results is None:
        raise RuntimeError("Could not find applet {}".format(first_applet))
    dxni_call_with_path(project, results["id"], version_id, verbose)


def native_call_app_setup(version_id, verbose):
    app_name = "native_hello"

    # Check if they already exist
    apps = list(dxpy.bindings.search.find_apps(name=app_name))
    if len(apps) == 0:
        # build the app
        cmdline = [
            "dx",
            "build",
            "--create-app",
            "--publish",
            os.path.join(top_dir, "test/apps/{}".format(app_name)),
        ]
        print(" ".join(cmdline))
        subprocess.check_output(cmdline)

    # build WDL wrapper tasks in test/dx_extern.wdl
    header_file = os.path.join(top_dir, "test/wdl_1_0/dx_app_extern.wdl")
    cmdline = [
        "java",
        "-jar",
        os.path.join(top_dir, "dxCompiler-{}.jar".format(version_id)),
        "dxni",
        "-apps",
        "only",
        "-force",
        "-language",
        "wdl_v1.0",
        "-output",
        header_file,
    ]
    if verbose:
        cmdline.append("--verbose")
    print(" ".join(cmdline))
    subprocess.check_output(cmdline)

    # check if providing an app-id in the path argument works
    results = dxpy.bindings.search.find_one_app(
        name=app_name, zero_ok=True, more_ok=False
    )
    if results is None:
        raise RuntimeError("Could not find app {}".format(app_name))
    dxni_call_with_path(None, results["id"], version_id, verbose)


######################################################################
# Compile the WDL files to dx:workflows and dx:applets
# delay_compile_errors: whether to aggregate all compilation errors
#   and only raise an Exception after trying to compile all the tests
def compile_tests_to_project(
    trg_proj,
    test_names,
    applet_folder,
    compiler_flags,
    version_id,
    lazy_flag,
    delay_compile_errors=False,
):
    runnable = {}
    has_errors = False
    for tname in test_names:
        oid = None
        if lazy_flag:
            oid = lookup_dataobj(tname, trg_proj, applet_folder)
        if oid is None:
            c_flags = compiler_flags[:] + compiler_per_test_flags(tname)
            try:
                oid = build_test(tname, trg_proj, applet_folder, version_id, c_flags)
            except subprocess.CalledProcessError:
                if delay_compile_errors:
                    traceback.print_exc()
                    has_errors = True
                else:
                    raise
        runnable[tname] = oid
        print("runnable({}) = {}".format(tname, oid))
    if has_errors:
        raise RuntimeError("failed to compile one or more tests")
    return runnable


def main():
    global test_unlocked
    argparser = argparse.ArgumentParser(
        description="Run WDL compiler tests on the platform"
    )
    argparser.add_argument(
        "--archive", help="Archive old applets", action="store_true", default=False
    )
    argparser.add_argument(
        "--compile-only",
        help="Only compile the workflows, don't run them",
        action="store_true",
        default=False,
    )
    argparser.add_argument("--compile-mode", help="Compilation mode")
    argparser.add_argument(
        "--debug",
        help="Run applets with debug-hold, and allow ssh",
        action="store_true",
        default=False,
    )
    argparser.add_argument(
        "--delay-workspace-destruction",
        help="Run applets with delayWorkspaceDestruction",
        action="store_true",
        default=False,
    )
    argparser.add_argument(
        "--force",
        help="Remove old versions of applets and workflows",
        action="store_true",
        default=False,
    )
    argparser.add_argument(
        "--folder", help="Use an existing folder, instead of building dxCompiler"
    )
    argparser.add_argument(
        "--lazy",
        help="Only compile workflows that are unbuilt",
        action="store_true",
        default=False,
    )
    argparser.add_argument(
        "--list",
        "--test-list",
        help="Print a list of available tests",
        action="store_true",
        dest="test_list",
        default=False,
    )
    argparser.add_argument(
        "--clean",
        help="Remove build directory in the project after running tests",
        action="store_true",
        default=False,
    )
    argparser.add_argument(
        "--delay-compile-errors",
        help="Compile all tests before failing on any errors",
        action="store_true",
        default=False,
    )
    argparser.add_argument(
        "--failed",
        help="Run the tests that failed previously (requires a .failed file in the current directory)",
        action="store_true",
        default=False
    )
    argparser.add_argument(
        "--locked",
        help="Generate locked-down workflows",
        action="store_true",
        default=False,
    )
    argparser.add_argument(
        "--project", help="DNAnexus project ID", default="dxCompiler_playground"
    )
    argparser.add_argument(
        "--project-wide-reuse",
        help="look for existing applets in the entire project",
        action="store_true",
        default=False,
    )
    argparser.add_argument(
        "--stream-all-files",
        help="Stream all input files with dxfs2",
        action="store_true",
        default=False,
    )
    argparser.add_argument(
        "--runtime-debug-level",
        help="printing verbosity of task/workflow runner, {0,1,2}",
    )
    argparser.add_argument(
        "--test", help="Run a test, or a subgroup of tests", action="append", default=[]
    )
    argparser.add_argument(
        "--unlocked",
        help="Generate only unlocked workflows",
        action="store_true",
        default=False,
    )
    argparser.add_argument(
        "--verbose", help="Verbose compilation", action="store_true", default=False
    )
    argparser.add_argument(
        "--verbose-key", help="Verbose compilation", action="append", default=[]
    )
    args = argparser.parse_args()

    print("top_dir={} test_dir={}".format(top_dir, test_dir))

    register_all_tests(args.verbose)
    if args.test_list:
        print_test_list()
        exit(0)
    test_names = []
    if args.failed and os.path.exists(".failed"):
        with open(".failed", "rt") as inp:
            test_names = [t.strip() for t in inp.readlines()]
    elif len(args.test) == 0:
        test_names = choose_tests("M")
    else:
        for t in args.test:
            test_names += choose_tests(t)
    print("Running tests {}".format(test_names))
    version_id = util.get_version_id(top_dir)

    project = util.get_project(args.project)
    if project is None:
        raise RuntimeError("Could not find project {}".format(args.project))
    if args.folder is None:
        base_folder = util.create_build_dirs(project, version_id)
    else:
        # Use existing prebuilt base folder
        base_folder = args.folder
        util.create_build_subdirs(project, base_folder)
    applet_folder = base_folder + "/applets"
    test_folder = base_folder + "/test"
    print("project: {} ({})".format(project.name, project.get_id()))
    print("folder: {}".format(base_folder))

    test_dict = {"aws:us-east-1": project.name + ":" + base_folder}

    # build the dxCompiler jar file, only on us-east-1
    assets = util.build(project, base_folder, version_id, top_dir, test_dict)
    print("assets: {}".format(assets))

    if args.unlocked:
        # Disable all locked workflows
        args.locked = False
        test_unlocked = test_names

    compiler_flags = []
    if args.locked:
        compiler_flags.append("-locked")
        test_unlocked = set()
    if args.archive:
        compiler_flags.append("-archive")
    if args.compile_mode:
        compiler_flags += ["-compileMode", args.compile_mode]
    if args.force:
        compiler_flags.append("-force")
    if args.verbose:
        compiler_flags.append("-verbose")
    if args.stream_all_files:
        compiler_flags.append("-streamAllFiles")
    if args.verbose_key:
        for key in args.verbose_key:
            compiler_flags += ["-verboseKey", key]
    if args.runtime_debug_level:
        compiler_flags += ["-runtimeDebugLevel", args.runtime_debug_level]
    if args.project_wide_reuse:
        compiler_flags.append("-projectWideReuse")

    #  is "native" included in one of the test names?
    if "call_native" in test_names or "call_native_v1" in test_names:
        native_call_setup(project, applet_folder, version_id, args.verbose)
    if "call_native_app" in test_names:
        native_call_app_setup(version_id, args.verbose)

    try:
        # Compile the WDL files to dx:workflows and dx:applets
        runnable = compile_tests_to_project(
            project,
            test_names,
            applet_folder,
            compiler_flags,
            version_id,
            args.lazy,
            args.delay_compile_errors,
        )
        if not args.compile_only:
            run_test_subset(
                project,
                runnable,
                test_folder,
                args.debug,
                args.delay_workspace_destruction,
            )
    finally:
        if args.clean:
            project.remove_folder(base_folder, recurse=True, force=True)
        print("Completed running tasks in {}".format(args.project))


if __name__ == "__main__":
    main()
