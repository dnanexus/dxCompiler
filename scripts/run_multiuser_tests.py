#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys
import time

import dxpy
from dxpy.exceptions import DXJobFailureError

import util

HERE = os.path.dirname(sys.argv[0])
TOP_DIR = os.path.dirname(os.path.abspath(HERE))
TEST_DIR = os.path.join(os.path.abspath(TOP_DIR), "test")

ALICE_SECURITY_CONTEXT = {
    "auth_token_type": "Bearer"
}
BOB_SECURITY_CONTEXT = {
    "auth_token_type": "Bearer"
}
VERSION_ID = ""
PROJECT = None
APPLET_FOLDER = ""
TEST_FOLDER = ""
BILLING_ORG = "org-dnanexus_apps"

def register_tokens(alice_token, bob_token):
    global ALICE_SECURITY_CONTEXT
    ALICE_SECURITY_CONTEXT["auth_token"] = alice_token
    global BOB_SECURITY_CONTEXT
    BOB_SECURITY_CONTEXT["auth_token"] = bob_token

def register_version(version_id):
    global VERSION_ID
    VERSION_ID = version_id

def register_project(project):
    global PROJECT
    PROJECT = project

def register_folders(applet_folder, test_folder):
    global APPLET_FOLDER
    APPLET_FOLDER = applet_folder
    global TEST_FOLDER
    TEST_FOLDER = test_folder

def login_alice():
    dxpy.set_api_server_info(host="stagingapi.dnanexus.com")
    dxpy.set_security_context(ALICE_SECURITY_CONTEXT)
    dxpy.set_project_context(PROJECT.get_id())
    print("Logged in as {}".format(dxpy.bindings.whoami()))

def login_bob():
    dxpy.set_api_server_info(host="stagingapi.dnanexus.com")
    dxpy.set_security_context(BOB_SECURITY_CONTEXT)
    dxpy.set_project_context(PROJECT.get_id())
    print("Logged in as {}".format(dxpy.bindings.whoami()))

def specific_applet_folder(tname):
    return "{}/{}".format(APPLET_FOLDER, tname)

def specific_test_folder(tname):
    return "{}/{}".format(TEST_FOLDER, tname)

def wait_for_completion(exec_objs):
    successes = []
    failures = []
    for i, exec_obj in exec_objs:
        dx_desc = exec_obj.describe()
        exec_name = dx_desc["name"].split(" ")[0]
        print("Awaiting completion {}.{}".format(exec_name, i))
        try:
            exec_obj.wait_on_done()
            print("Analysis succeeded {}.{}".format(exec_name, i))
            successes.append((exec_name, exec_obj, i))
        except DXJobFailureError:
            print("Error: analysis failed {}.{}".format(exec_name, i))
            failures.append((exec_name, exec_obj, i))

    return successes, failures

def report_test_success(tname):
    print("Test {} succeeded".format(tname))

def report_test_failure(tname):
    print("Test {} failed".format(tname))

def test_global_wf_from_wdl():
    # As Alice, compile workflow and create global workflow
    login_alice()

    tname = "global_wf_from_wdl"
    wf_source = os.path.join(os.path.abspath(TEST_DIR), "multi_user", "{}.wdl".format(tname))
    compiler_flags = ["-instanceTypeSelection", "dynamic"]
    workflow_id = util.build_executable(
        source_file=wf_source,
        project=PROJECT,
        folder=specific_applet_folder(tname),
        top_dir=TOP_DIR,
        version_id=VERSION_ID,
        compiler_flags=compiler_flags)
    full_workflow_id = "{}:{}".format(PROJECT.get_id(), workflow_id)
    print("Compiled {}".format(full_workflow_id))

    # Strictly increasing global WF version
    global_workflow_version = "1.0.{}".format(int(time.time()))
    print("Global workflow version {}".format(global_workflow_version))

    # Build global workflow from workflow
    global_workflow_name = "globalworkflow-{}".format(tname)
    print("Global workflow name {}".format(global_workflow_name))
    global_wf_build_cmd = [
        "dx",
        "build",
        "--globalworkflow",
        "--from",
        full_workflow_id,
        "--version",
        global_workflow_version,
        "--bill-to",
        BILLING_ORG
    ]
    global_wf_publish_cmd = [
        "dx",
        "publish",
        "{}/{}".format(global_workflow_name, global_workflow_version)
    ]
    try:
        print(" ".join(global_wf_build_cmd))
        subprocess.check_call(global_wf_build_cmd)
        print(" ".join(global_wf_publish_cmd))
        subprocess.check_call(global_wf_publish_cmd)
    except subprocess.CalledProcessError as cpe:
        print("Error building global workflow from {}\n stdout: {}\n stderr: {}".format(
            workflow_id,
            cpe.stdout,
            cpe.stderr
        ))
        raise

    # Do some developer actions on global workflow
    add_developers_cmd = [
        "dx",
        "add",
        "developers",
        global_workflow_name,
        "org-dnanexus_apps"
    ]
    add_users_cmd = [
        "dx",
        "add",
        "users",
        global_workflow_name,
        "user-dnanexus_apps_test_robot"
    ]
    add_categories_cmd = [
        "dx",
        "api",
        global_workflow_name,
        "addCategories",
        '{"categories":["category_1", "category_2"]}'
    ]
    remove_categories_cmd = [
        "dx",
        "api",
        global_workflow_name,
        "removeCategories",
        '{"categories":["category_1", "category_2"]}'
    ]
    add_tags_cmd = [
        "dx",
        "api",
        "{}/{}".format(global_workflow_name, global_workflow_version),
        "addTags",
        '{"tags":["tag_1", "tag_2"]}'
    ]
    remove_tags_cmd = [
        "dx",
        "api",
        "{}/{}".format(global_workflow_name, global_workflow_version),
        "removeTags",
        '{"tags":["tag_1", "tag_2"]}'
    ]
    update_metadata_cmd = [
        "dx",
        "api",
        "{}/{}".format(global_workflow_name, global_workflow_version),
        "update",
        '{"title":"Global WF from WDL 1", "summary":"Summary 1", "developerNotes":"Notes 1"}'
    ]

    try:
        print(" ".join(add_developers_cmd))
        subprocess.check_call(add_developers_cmd)
        print(" ".join(add_users_cmd))
        subprocess.check_call(add_users_cmd)
        print(" ".join(add_categories_cmd))
        subprocess.check_call(add_categories_cmd)
        print(" ".join(remove_categories_cmd))
        subprocess.check_call(remove_categories_cmd)
        print(" ".join(add_tags_cmd))
        subprocess.check_call(add_tags_cmd)
        print(" ".join(remove_tags_cmd))
        subprocess.check_call(remove_tags_cmd)
        print(" ".join(update_metadata_cmd))
        subprocess.check_call(update_metadata_cmd)
    except subprocess.CalledProcessError as cpe:
        print("Error during developer actions on {}\n stdout: {}\n stderr: {}".format(
            global_workflow_name,
            cpe.stdout,
            cpe.stderr
        ))
        raise


"""Commented out in APPS-1197 because of the fund limits on Bob's account
    # As Bob, run global workflow
    login_bob()

    exec_objs = util.run_executable(
        oid=global_workflow_name,
        project=PROJECT,
        test_folder=specific_test_folder(tname),
        test_name=tname
    )
"""
    successes, failures = wait_for_completion(exec_objs)
    if len(successes) > 0:
        report_test_success(tname)
    elif len(failures) > 0:
        report_test_failure(tname)
        raise RuntimeError("Analysis failed in test {}".format(tname))

def main():
    argparser = argparse.ArgumentParser(
        description="Run dxCompiler multi-user tests on the platform"
    )
    argparser.add_argument(
        "--alice-token",
        help="Token for user-dnanexus_apps_robot on staging",
        required=True,
    )
    argparser.add_argument(
        "--bob-token",
        help="Token for user-dnanexus_apps_test_robot on staging",
        required=True,
    )
    argparser.add_argument(
        "--folder", help="Use an existing folder with dxCompiler assets, instead of building dxCompiler"
    )
    argparser.add_argument(
        "--project", help="DNAnexus project ID", default="dxCompiler_playground"
    )
    argparser.add_argument(
        "--clean",
        help="Remove build directory in the project after running tests",
        action="store_true",
        default=False,
    )
    args = argparser.parse_args()

    register_tokens(args.alice_token, args.bob_token)

    version_id = util.get_version_id(TOP_DIR)
    register_version(version_id)

    project = util.get_project(args.project)
    if project is None:
        raise RuntimeError("Could not find project {}".format(args.project))
    print("Project {} ({})".format(project.name, project.get_id()))
    register_project(project)

    # Do folder setup as Alice
    login_alice()

    if args.folder is None:
        base_folder = util.create_build_dirs(project, version_id)
    else:
        # Use an existing folder with dxCompiler assets, instead of building dxCompiler
        base_folder = args.folder
        util.create_build_subdirs(project, base_folder)
    print("Base folder {}".format(base_folder))
    applet_folder = base_folder + "/applets"
    test_folder = base_folder + "/test"
    register_folders(applet_folder, test_folder)

    # Build the dxCompiler jar file, only on us-east-1
    test_dict = {"aws:us-east-1": project.name + ":" + base_folder}
    assets = util.build(project, base_folder, version_id, TOP_DIR, test_dict,
                        force=False)
    print("assets: {}".format(assets))

    # TODO APPS-1030 Improve this script by using pytest

    # Run tests
    try:
        test_global_wf_from_wdl()
    finally:
        if args.clean:
            # Do folder cleanup as Alice
            login_alice()
            project.remove_folder(base_folder, recurse=True, force=True)
        print("Completed running tasks in {}".format(args.project))

if __name__ == "__main__":
    main()
