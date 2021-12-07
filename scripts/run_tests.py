#!/usr/bin/env python3
import argparse
import shutil
from collections import namedtuple
import dxpy
import glob
import hashlib
import json
import os
import random
import re
import sys
import subprocess
import tempfile
from termcolor import cprint
import time
import traceback
from typing import List
import yaml
from dxpy.exceptions import DXJobFailureError

import util

here = os.path.dirname(sys.argv[0])
top_dir = os.path.dirname(os.path.abspath(here))
test_dir = os.path.join(os.path.abspath(top_dir), "test")
default_instance_type = "mem1_ssd1_v2_x4"

git_revision = subprocess.check_output(
    ["git", "describe", "--always", "--dirty", "--tags"]
).strip()
test_files = {}

expected_failure = {
    "bad_status",
    "bad_status2",
    "just_fail_wf",
    "missing_output",
    "docker_retry",
    "argument_list_too_long",
    "diskspace_exhauster",
    "symlink-illegal",
    "docker-array-secondaryfiles.1",
    "iwd-container-entryname2",
    "iwd-container-entryname3",
    "iwd-container-entryname4",
    "loadContents-limit",
    "cond-wf-003.3",
    "cond-wf-004.1",
    "cond-wf-005",
    "cond-wf-006.1",
    "cond-wf-012",
    "cond-wf-003-1.1",
    "cond-wf-003-1_nojs.1",
    "cond-wf-004_nojs.1",
    "cond-wf-005_nojs",
    "cond-wf-006_nojs.1",
    "cond-wf-012_nojs",
    "fail-unconnected",
}

test_compilation_failing = {"import_passwd"}

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
    "top_wf",
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
    "apps_936",

    # manifests
    "simple_manifest",
    "complex_manifest",
    "view_and_count_manifest",

    # workflow with output files created by expressions
    "upload_workflow_files",
    "subworkflow_with_task",
    "apps_700",
    "apps_864"
]

wdl_v1_1_list = [
    "v1_1_dict",
    "apps_847_scatter_empty",
    "optional_missing",
    "inputs_provided_optional",

    # bug regression tests
    "apps_579_boolean_flag_expr",
    "apps_579_string_substitution_expr",
    "apps_956_private_var_local"
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
    # "population",

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
    os.path.basename(path)[:-9]
    for path in glob.glob(os.path.join(test_dir, "cwl_conformance", "tools", "*.cwl.json"))
]
cwl_conformance_workflows = [
    os.path.basename(path)[:-9]
    for path in glob.glob(os.path.join(test_dir, "cwl_conformance", "workflows", "*.cwl.json"))
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

cromwell_key_error_list = [
    "http_inputs",
    "drs_usa_hca",
    "drs_usa_jdr",
]

# These are cromwell tests that won't run on DNAnexus - see README.txt
cromwell_invalid = {
    "local_backend",
    "string_interpolation",
    "call_cache_hit_prefixes",
    "declarations",
    "reference_disk_test",
    "optional_parameter",
    "sub",
    "sub_sub",
    "echo",
    "sub_workflow_no_output",
    "recursive_imports",
    "large_final_workflow_outputs_dir",
    "input_from_bucket_with_requester_pays",
    "optional_declarations",
    "sub_workflow_interactions",
    "unscattered",
    "inter_scatter_dependencies",
    "docker_alpine",
    "read_write_functions",
    "afters_and_ifs",
    "afters",
    "afters_and_scatters",
    "custom_cacheworthy_attributes",
    "input_expressions",
    "missing_delete",
    "confirm_preemptible",
    "call_cache_capoeira_jes",
    "dedup_localizations_papi_v2",
    "papi_v2_log",
    "papi_v2_plain_detritus",
    "call_cache_capoeira_local",
    "backendWithNoDocker",
    "docker_image_cache_true",
    "dummy_scatter",
    "fofn_caching",
    "hello_private_repo",
    "local_bourne",
    "papi_v2_gcsa",
    "monitoring_log",
    "call_cache_capoeira_tes",
    "check_network_in_vpc",
    "tmp_dir",
    "long_cmd",
    "workbench_health_monitor_check",
    "monitoring_image_script",
    "docker_size_dockerhub",
    "docker_size_gcr",
    "custom_mount_point",
    "short_circuit",
    "top",
    "recursive_imports_no_subwf",
    "parallel_composite_uploads_on",
    "parallel_composite_uploads_off",
    "default_runtime_attributes",
}

# tests taken from cromwell repository
cromwell_tests_list = [
    "null_input_values",
    "dont_strip_line_prefix",
    "non_root_default_user",
    "memory_units",
    "cacheWithinWF",
    "dot_dir_stuck_running",
    "empty_string",
    "floating_tags",
    "array_literal_locations",
    "stdout_delete",
    "sub_workflow_delete",
    "no_output_delete",
    "exhaustive_delete",
    "scatter_delete",
    "collections_delete",
    "hello_delete",
    "sub_workflow_delete_import",
    "no_cache_delete",
    "gcs_path_ending_with_newline",
    "readFromCache",
    "sizerelativepath",
    "subworkflow_wt",
    "b",
    "c",
    "a",
    "d",
    "sub_sub_sub",
    "array_io",
    "simple_if",
    "single_to_array_conversion",
    "coerce_to_array_of_optionals",
    "wdl_function_locations",
    "workflow_output_paths",
    "sub_function",
    "public_http_import",
    "control_chars",
    "prefix",
    "write_lines_files",
    "cached_copy",
    "read_tsv",
    "custom_entrypoint",
    "square",
    "papi_cpu_platform",
    "complex_types_files",
    "file_evaluator_identifier_lookups",
    "non_root_specified_user",
    "write_lines",
    "workflow_output_paths_colliding",
    "jes_labels",
    "localization_sanity_papi_v2",
    "recimp_nosubwf_outer",
    "recimp_nosubwf_inner",
    "globbingindex",
    "postfix_quantifiers",
    "length",
    "wdl_empty_glob",
    "output_filename_interpolation",
    "aliased_subworkflows",
    "docker_image_cache_false",
    "curl",
    "symlink_localization",
    "error_10_preemptible",
    "multiline_command_line",
    "use_cacheCopy_dir",
    "writeToCache",
    "cacheBetweenWF",
    "lots_of_inputs",
    "local_gcs",
    "read_write_json_roundtrip_develop",
    "read_write_json_roundtrip",
    "checkpointing",
    "cromwell_restart",
    "space",
    "arrays_scatters_ifs",
    "declarations_as_nodes",
    "variable_scoping",
    "sub_workflow_decls",
    "input_mirror",
    "sub_workflow_hello_world_import",
    "sub_workflow_hello_world",
    "volatile_disables_cache",
    "file_outputs_from_input",
    "write_tsv",
    "final_call_logs_dir",
    "subdirectory",
    "input_localization",
    "scattered",
    "filearrayoutput",
    "array_io",
    "docker_hash_quay",
    "docker_hash_gcr",
    "workflow_type_and_version_wdl",
    "dontglobinputs",
    "globbingscatter",
    "ifs_in_scatters",
    "nested_lookups",
    "simple_if",
    "declarations_in_ifs",
    "lots_of_nesting",
    "ifs_upstream_and_downstream",
    "subworkflows_in_ifs",
    "scatters_in_ifs",
    "simple_if_workflow_outputs",
    "scattergather",
    "map_workflow",
    "forkjoin",
    "scatter_chain",
    "output_redirection",
    "workflowenginefunctions",
    "stdout_stderr_passing",
    "scatter",
    "siblings_scatter",
    "simple_scatter",
    "prepare_scatter_gather",
    "multiplesourcedarray",
    "passingfiles",
    "referencingpreviousinputsandoutputs",
    "engine_functions",
    #"string_interpolation_optional",  # pending wdlTools 170
    #"none_literal",  # pending wdlTools 170
    "sub_workflow_interactions_scatter",
    "sub_workflow_one_output_import",
    "sub_workflow_var_refs",
    "sub_workflow_var_refs_import",
    #"globbingBehavior",  # pending dxCompiler 87
    #"object_access",  # pending wdlTools 171
    #"read_write_json",  # pending wdlTools 171
    "no_task_no_output_delete",
    "if_then_else_expressions",
    "sub_workflow_no_output_block_import",
    "sub_workflow_no_outputs_in_block_import",
    "sub_workflow_interactions_import",
    "workflow_output_declarations",
    "member_access",
    "select_functions",
    "dollars_in_strings",
    "workflow_name_length_ok",
    "importer_ok",
    "read_write_map",
    "docker_image_cache_unspecified",
    "defined_function",
    "workflow_engine_functions",
    "empty_scatter",
    "continue_on_return_code",
    "exit",
]

cwl_cromwell_tests_list = [
    "cwl_ad_hoc_file_test",
    "cwl_cache_between_workflows",
    "cwl_cache_within_workflow",
    "cwl_docker_size",
    "cwl_dynamic_initial_workdir",
    "cwl_expressionLib",
    "cwl_format",
    # "cwl_format_url", # Could not load extension schema https
    "cwl_glob_sort",
    "cwl_hello",
    # "cwl_http_inputs", # HTTPS input link is not supported: 
    #                     Error translating inputs: java.lang.RuntimeException: Unsupported file source .png
    "test_wf",
    "touch",
    "test_pack",
    "cwl_input_binding_expression",
    # "cwl_input_json", # Error translating to IR
    "cwl_input_typearray",
    "cwl_interpolated_strings",
    "cwl_optionals",
    "cwl_output_json",
    "prefix_for_array",
    "cwl_recursive_link_directories",
    "cwl_relative_imports",
    # "cwl_disk_resources", # Could not resolve host: metadata.google.internal
    #                         # Unknown hint https://www.dnanexus.com/cwl#InputResourceRequirement (Should be deprecated)
    # "cwl_inputdir_zero_doesnt_localize", # Error translating to IR
    # "cwl_resources", # Could not resolve host: metadata.google.internal
    # "cwl_restart", # AppInternalError: workflow does not contain a tool APPS-834
    "1st-tool",
    "cwl_secondary_files",
    # "cwl_secondary_files_workflow", # Error creating translator
    "cwl_stdout_expression",
    # "scatter-wf1", # Could not find linking information APPS-834
    # "cwl_three_step", # AppInternalError: workflow does not contain a tool APPS-834
    # "cwl_three_step_caller_wf" # Error translating to IR
]

# these are tests that take a long time to run
long_test_list = [
    "diskspace_exhauster"  # APPS-749
]

medium_test_list = (
    wdl_v1_list + wdl_v1_1_list + docker_test_list + special_flags_list + cwl_tools
)
large_test_list = (
    medium_test_list + draft2_test_list + single_tasks_list + doc_tests_list + long_test_list +
    cwl_conformance_tools + cwl_conformance_workflows + cromwell_tests_list
)

manifest_test_list = ("simple_manifest", "complex_manifest", "view_and_count_manifest")

test_suites = {
    "CI": ci_test_list,
    "M": medium_test_list,
    "L": large_test_list,
    "tasks": single_tasks_list,
    "draft2": draft2_test_list,
    "docker": docker_test_list,
    "native": ["call_native", "call_native_v1"],
    "docs": doc_tests_list,
    "cwl_tools": cwl_conformance_tools,
    "cwl_workflows": cwl_conformance_workflows,
    'cromwell': cromwell_tests_list,
    "cwl_cromwell": cwl_cromwell_tests_list,
    'manifests': manifest_test_list
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

# use the applet's default instance type rather than the default (mem1_ssd1_x4)
test_instance_type = ["diskspace_exhauster"]


# Read a JSON file
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
        # the id in a packed CWL file is always "main" so we use the test name instead
        if doc["class"] == "Workflow":
            return TestMetaData(name=tname, kind="workflow")
        elif doc["class"] == "CommandLineTool":
            return TestMetaData(name=tname, kind="applet")
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
    elif ext == ".cwl.json":
        metadata = get_cwl_json_metadata(source_file, tname)
    elif ext == ".cwl":
         metadata = get_cwl_metadata(source_file, tname)
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
    elif os.path.exists(os.path.join(dir_path, tname + "_input.yaml")):
        test_yaml = os.path.join(dir_path, tname + "_input.yaml")
        desc.raw_input.append(test_yaml)
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
        return dxpy.DXFile(fields, project.id)
    else:
        return dxpy.DXFile(fields["id"], fields.get("project", project.id))


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
        expected_checksum = expected_val["checksum"]
        algo = expected_checksum.split("$")[0]
    else:
        algo = None
        expected_checksum = None

    location = None
    size = None
    checksum = None
    secondary_files = None

    if isinstance(result, str):
        contents = result.strip()
    elif result.get("class", result.get("type")) == "File":
        contents = result.get("contents")
        location = result.get("location", result.get("path", result.get("uri")))
        if isinstance(location, dict) and "$dnanexus_link" in location:
            dxfile = link_to_dxfile(location, project)
            location = os.path.join(dxfile.describe()["folder"], dxfile.describe()["name"])
            if contents is None:
                contents = download_dxfile(dxfile)
        size = result.get("size")
        checksum = result.get("checksum")
        secondary_files = result.get("secondaryFiles", [])
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
        if size is None and contents is not None:
            size = len(contents)
        if size != expected_val["size"]:
            if verbose:
                cprint("Analysis {} gave unexpected results".format(tname), "red")
                cprint(
                    "Field {} should have size ({}), actual: ({})".format(
                        field_name, expected_val["size"], size
                    ),
                    "red"
                )
            return False

    if expected_checksum:
        checksum = checksum or (get_checksum(contents, algo) if contents else None)
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
        if len(expected_secondary_files) != seconary_files_len:
            if verbose:
                cprint("Analysis {} gave unexpected results".format(tname), "red")
                cprint(
                    "Field {} should have secondaryFiles ({}), actual = ({})".format(
                        field_name, expected_secondary_files, secondary_files
                    ),
                    "red"
                )
            return False
        # TODO: sort both lists rather than doing an all-by-all comparison
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


# expected_val is a serialized CwlDirectory (an object with class="Directory")
# result may be a serialized CwlDirectory or a serialized VFolder (an object
#   with type="Folder")
def compare_result_directory(result, expected_val, field_name, tname, project, verbose=True):
    location = result.get("location", result.get("path", result.get("uri")))
    if location is not None and location.startswith("dx://"):
        project_id, folder = location[5:].split(":")
        project = dxpy.DXProject(project_id)
    else:
        folder = location

    if "basename" in result:
        basename = result["basename"]
    elif folder:
        basename = os.path.basename(folder)
    else:
        basename = None

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
    expected_cls = expected_val.get("class", expected_val.get("type"))
    if cls == "File" and expected_cls == "File":
        return compare_result_file(result, expected_val, field_name, tname, project, verbose)
    elif cls in {"Directory", "Folder"} and expected_cls in {"Directory", "Folder"}:
        return compare_result_directory(result, expected_val, field_name, tname, project, verbose)
    else:
        cprint("Analysis {} gave unexpected results".format(tname), "red")
        cprint(
            "Field {} should be of class ({}), actual = ({})".format(
                field_name, expected_cls, cls
            ),
            "red",
        )
        return False


# Check that a workflow returned the expected result for
# a [key]
def validate_result(tname, exec_outputs: dict, key, expected_val, project):
    if exec_outputs is None:
        if expected_val is None:
            return True
        else:
            cprint("Outputs missing for {}".format(tname), "red")
            return False

    desc = test_files[tname]
    # Extract the key. For example, for workflow "math" returning
    # output "count":
    #    'math.count' -> count
    (exec_name, *field_name_parts) = key.split(".")

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
        elif expected_val is None:
            # optional
            return True
        else:
            cprint(
                "Field {} missing from executable results {}".format(
                    field_name1, exec_outputs
                ),
                "red",
            )
            return False

        # Sort two lists of dicts to make them comparable. Given lists of dicts a and b:
        # 1. get the set of all keys in all dicts in both lists
        # 2. expand each dict into a list of tuples where the first element is the key and the
        # second element is the value (which may be None if the dict does not contain the key)
        # 3. sort each list
        # 4. remove the tuples with None values from each list
        # 5. turn each list of tuples back into a list of dicts.
        def sort_dicts(a, b):
            all_keys = list(sorted(
                set(k for x in a for k in x.keys()) | set(k for x in b for k in x.keys())
            ))
            return (
                [
                    dict((k, v) for k, v in x if v is not None)
                    for x in list(sorted([(k, i.get(k)) for k in all_keys] for i in a))
                ],
                [
                    dict((k, v) for k, v in x if v is not None)
                    for x in list(sorted([(k, j.get(k)) for k in all_keys] for j in b))
                ]
            )

        def sort_maybe_mixed(seq):
            try:
                # may fail if the lists contain mutliple types of values
                return list(sorted(seq))
            except:
                d = dict((str(x), x) for x in seq)
                sorted_keys = list(sorted(d.keys()))
                return [d[k] for k in sorted_keys]

        def compare_values(expected, actual, field):
            if isinstance(actual, dict) and "___" in actual:
                actual = actual["___"]
                if isinstance(expected, dict) and "___" in expected:
                    expected = expected["___"]
            if isinstance(actual, dict) and "wrapped___" in actual:
                actual = actual["wrapped___"]
                if isinstance(expected, dict) and "wrapped___" in expected:
                    expected = expected["wrapped___"]

            if isinstance(actual, list) and isinstance(expected, list):
                actual = list(filter(lambda x: x is not None, actual))
                expected = list(filter(lambda x: x is not None, expected))
                n = len(actual)
                if n != len(expected):
                    cprint("Analysis {} gave unexpected results".format(tname), "red")
                    cprint(
                        "Field {} should have length ({}), actual = ({})".format(
                            field, len(expected), len(actual)
                        ),
                        "red",
                    )
                    return False
                if n == 0:
                    return True
                elif n > 1:
                    if isinstance(actual[0], dict):
                        actual, expected = sort_dicts(actual, expected)
                    else:
                        actual = sort_maybe_mixed(actual)
                        expected = sort_maybe_mixed(expected)

                for i, (e, a) in enumerate(zip(expected, actual)):
                    if not compare_values(e, a, "{}[{}]".format(field, i)):
                        return False
                else:
                    return True

            if isinstance(expected, dict) and (
                expected.get("class") in {"File", "Directory"} or
                expected.get("type") in {"File", "Folder"}
            ):
                return compare_result_path(actual, expected, field_name1, tname, project)

            if isinstance(actual, dict) and actual.get("type") == "File" and "uri" in actual:
                actual = actual["uri"]
            if isinstance(actual, dict) and "$dnanexus_link" in actual:
                actual = download_dxfile(link_to_dxfile(actual, project))

            if isinstance(actual, dict) and isinstance(expected, dict):
                expected_keys = set(expected.keys())
                actual_keys = set(expected.keys())
                if expected_keys != actual_keys:
                    cprint("Analysis {} gave unexpected results".format(tname), "red")
                    cprint(
                        "Field {} should have keys ({}), actual = ({})".format(
                            field, expected_keys, actual_keys
                        ),
                        "red",
                    )
                    return False
                for k in expected_keys:
                    if not compare_values(expected[k], actual[k], "{}[{}]".format(field, k)):
                        return False
                else:
                    return True

            if str(actual).strip() != str(expected).strip():
                cprint("Analysis {} gave unexpected results".format(tname), "red")
                cprint(
                    "Field {} should be ({}), actual = ({})".format(
                        field, expected, actual
                    ),
                    "red",
                )
                return False

            return True

        return compare_values(expected_val, result, field_name1)
    except:
        traceback.print_exc()
        return False


def get_checksum(contents, algo):
    try:
        m = hashlib.new(algo)
        m.update(contents)
        checksum = m.digest()
        return f"{algo}${checksum}"
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
    # both static and dynamic instance type selection should work,
    # so we can test them at random
    instance_type_selection = random.choice(["static", "dynamic"])
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
        "-instanceTypeSelection",
        instance_type_selection
    ]
    if "manifest" in desc.source_file:
        cmdline.append("-useManifests")
    cmdline += compiler_flags
    print(" ".join(cmdline))
    try:
        oid = subprocess.check_output(cmdline).strip()
    except subprocess.CalledProcessError as cpe:
        print(f"error compiling {desc.source_file}\n  stdout: {cpe.stdout}\n  stderr: {cpe.stderr}")
        raise
    return oid.decode("ascii")


def ensure_dir(path):
    print("making sure that {} exists".format(path))
    if not os.path.exists(path):
        os.makedirs(path)


def wait_for_completion(test_exec_objs):
    print("awaiting completion ...")
    successes = []
    failures = []
    for i, exec_obj in test_exec_objs:
        tname = find_test_from_exec(exec_obj)
        desc = test_files[tname]
        try:
            exec_obj.wait_on_done()
            print("Analysis {}.{} succeeded".format(desc.name, i))
            successes.append((i, exec_obj, True))
        except DXJobFailureError:
            if tname in expected_failure or "{}.{}".format(tname, i) in expected_failure:
                print("Analysis {}.{} failed as expected".format(desc.name, i))
                successes.append((i, exec_obj, False))
            else:
                cprint("Error: analysis {}.{} failed".format(desc.name, i), "red")
                failures.append((tname, exec_obj))
    print("tools execution completed")
    return successes, failures


# Run [workflow] on several inputs, return the analysis ID.
def run_executable(
    project, test_folder, tname, oid, debug_flag, delay_workspace_destruction, instance_type=default_instance_type
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
                run_kwargs = {"ignore_reuse_stages": ["*"]}
            elif desc.kind == "applet":
                exec_obj = dxpy.DXApplet(project=project.get_id(), dxid=oid)
                run_kwargs = {"ignore_reuse": True}
            else:
                raise RuntimeError("Unknown kind {}".format(desc.kind))

            if debug_flag:
                run_kwargs["debug"] = {
                    "debugOn": ["AppError", "AppInternalError", "ExecutionError"]
                }
                run_kwargs["allow_ssh"] = ["*"]

            if delay_workspace_destruction:
                run_kwargs["delay_workspace_destruction"] = True
            if instance_type:
                run_kwargs["instance_type"] = instance_type

            return exec_obj.run(
                inputs,
                project=project.get_id(),
                folder=test_folder,
                name="{} {}".format(desc.name, git_revision),
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
        return [(0, run(-1))]
    else:
        return [(i, run(i)) for i in range(n)]


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
    project, runnable, test_folder, debug_flag, delay_workspace_destruction, delay_run_errors
):
    # Run the workflows
    test_exec_objs = []
    errors = [] if delay_run_errors else None
    for tname, oid in runnable.items():
        desc = test_files[tname]
        print("Running {} {} {}".format(desc.kind, desc.name, oid))
        if tname in test_instance_type:
            instance_type = None
        else:
            instance_type = default_instance_type
        try:
            anl = run_executable(
                project, test_folder, tname, oid, debug_flag, delay_workspace_destruction, instance_type
            )
            test_exec_objs.extend(anl)
        except Exception as ex:
            if tname in test_compilation_failing:
                cprint("Workflow {} execution failed as expected".format(tname))
                continue
            elif delay_run_errors:
                cprint(f"Workflow {tname} execution failed", "red")
                traceback.print_exc()
                errors.append(tname)
            else:
                raise ex

    if errors:
        write_failed(errors)
        raise RuntimeError(f"failed to run one or more tests {','.join(errors)}")

    print("executions: " + ", ".join([a[1].get_id() for a in test_exec_objs]))

    # Wait for completion
    successful_executions, failed_executions = wait_for_completion(test_exec_objs)

    print("Verifying results")

    def verify_test(exec_obj, i):
        exec_desc = exec_obj.describe()
        tname = find_test_from_exec(exec_obj)
        test_desc = test_files[tname]
        try:
            exec_outputs = extract_outputs(tname, exec_desc)
        except:
            if tname in expected_failure or "{}.{}".format(tname, i) in expected_failure:
                print("Analysis {}.{} failed as expected".format(tname, i))
                return None
            else:
                raise
        if len(test_desc.results) > i:
            shouldbe = read_json_file_maybe_empty(test_desc.results[i])
            correct = True
            print("Checking results for workflow {} job {}".format(test_desc.name, i))
            for key, expected_val in shouldbe.items():
                if not validate_result(tname, exec_outputs, key, expected_val, project):
                    correct = False
                    break
            if correct:
                if tname in expected_failure or "{}.{}".format(tname, i) in expected_failure:
                    cprint(
                        "Error: analysis {}.{} was expected to fail but its results are valid".format(test_desc.name, i),
                        "red"
                    )
                    return tname
                else:
                    print("Analysis {}.{} results are valid".format(test_desc.name, i))
                    return None
            else:
                if tname in expected_failure or "{}.{}".format(tname, i) in expected_failure:
                    print("Analysis {}.{} results are invalid as expected".format(test_desc.name, i))
                    return None
                else:
                    cprint("Error: analysis {}.{} results are invalid".format(test_desc.name, i), "red")
                    return tname

    failed_verifications = []
    for i, exec_obj, verify in successful_executions:
        if verify:
            failed_name = verify_test(exec_obj, i)
            if failed_name is not None:
                failed_verifications.append(failed_name)

    print("-----------------------------")
    print(f"Total tests: {len(test_exec_objs)}")

    if failed_executions or failed_verifications:
        failed_tools = set(e[0] for e in failed_executions)
        unverified_tools = set(failed_verifications)
        if failed_executions:
            fexec = "\n".join(failed_tools)
            print(f"Failed executions: {len(failed_executions)}")
            print(f"Tools failed execution:\n{fexec}")
        if failed_verifications:
            fveri = "\n".join(unverified_tools)
            print(f"Failed verifications: {len(failed_verifications)}")
            print(f"Tools failed results verification:\n{fveri}")
        write_failed(failed_tools | unverified_tools)
        raise RuntimeError("Failed")
    else:
        print("All tests successful!")


def write_failed(failed):
    # write failed tests to a file so we can easily re-run them next time
    # if a .failed file already exists, make a backup
    if os.path.exists(".failed"):
        bak_file = ".failed.bak"
        i = 0
        while os.path.exists(bak_file):
            bak_file = f".failed.bak.{i}"
            i += 1
        shutil.copy(".failed", bak_file)
    with open(".failed", "wt") as out:
        failed_sorted = sorted(set(tname.split(".")[0] for tname in failed))
        out.write("\n".join(failed_sorted))


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
    errors = [] if delay_compile_errors else None
    for tname in test_names:
        specific_applet_folder = "{}/{}".format(applet_folder, tname)
        oid = None
        if lazy_flag:
            oid = lookup_dataobj(tname, trg_proj, specific_applet_folder)
        if oid is None:
            c_flags = compiler_flags[:] + compiler_per_test_flags(tname)
            try:
                oid = build_test(tname, trg_proj, specific_applet_folder, version_id, c_flags)
            except subprocess.CalledProcessError:
                if tname in test_compilation_failing:
                    print("Workflow {} compilation failed as expected".format(tname))
                    continue
                elif delay_compile_errors:
                    cprint(f"Workflow {tname} compilation failed", "red")
                    traceback.print_exc()
                    errors.append(tname)
                else:
                    raise
        runnable[tname] = oid
        print("runnable({}) = {}".format(tname, oid))
    if errors:
        write_failed(errors)
        raise RuntimeError(f"failed to compile one or more tests: {','.join(errors)}")
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
        "--build",
        help="force: remove existing dxCompiler JAR and rebuild; only: only build dxCompiler, "
             "do not run any tests; if not specified, dxCompiler will be built only if there is "
             "not already a dxCompiler asset in the project",
        default=None
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
        "--delay-run-errors",
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
    elif args.test:
        for t in args.test:
            test_names += choose_tests(t)
    elif args.build != "only":
        test_names = choose_tests("M")
    if test_names:
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
    assets = util.build(project, base_folder, version_id, top_dir, test_dict,
                        force=args.build is not None)
    print("assets: {}".format(assets))

    if args.build == "only":
        exit(0)

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
                args.delay_run_errors,
            )
    finally:
        if args.clean:
            project.remove_folder(base_folder, recurse=True, force=True)
        print("Completed running tasks in {}".format(args.project))


if __name__ == "__main__":
    main()
