#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys
import time

import dxpy

import util

here = os.path.dirname(sys.argv[0])
top_dir = os.path.dirname(os.path.abspath(here))
test_dir = os.path.join(os.path.abspath(top_dir), "test")

ALICE_SECURITY_CONTEXT = {
    "auth_token_type": "Bearer"
}
BOB_SECURITY_CONTEXT = {
    "auth_token_type": "Bearer"
}
VERSION_ID = ""
PROJECT_ID = ""
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

def register_project(project_id):
    global PROJECT_ID
    PROJECT_ID = project_id

def register_folders(applet_folder, test_folder):
    global APPLET_FOLDER
    APPLET_FOLDER = applet_folder
    global TEST_FOLDER
    TEST_FOLDER = test_folder

def login_alice():
    dxpy.set_api_server_info(host="stagingapi.dnanexus.com")
    dxpy.set_security_context(ALICE_SECURITY_CONTEXT)
    dxpy.set_project_context(PROJECT_ID)
    print("Logged in as {}".format(dxpy.bindings.whoami()))

def login_bob():
    dxpy.set_api_server_info(host="stagingapi.dnanexus.com")
    dxpy.set_security_context(BOB_SECURITY_CONTEXT)
    dxpy.set_project_context(PROJECT_ID)
    print("Logged in as {}".format(dxpy.bindings.whoami()))

def specific_applet_folder(tname):
    return "{}/{}".format(APPLET_FOLDER, tname)

# Build a workflow.
#
# wf_source         Path to workflow source file
# project_id        Destination project id
# folder            Destination folder
# version_id        dxCompiler version
# compiler_flags    dxCompiler flags
def build_workflow(wf_source, project_id, folder, version_id, compiler_flags):
    print("Compiling {}".format(wf_source))
    cmdline = [
        "java",
        "-jar",
        os.path.join(top_dir, "dxCompiler-{}.jar".format(version_id)),
        "compile",
        wf_source,
        "-force",
        "-folder",
        folder,
        "-project",
        project_id
    ]
    cmdline += compiler_flags
    # print(" ".join(cmdline))
    # debug
    print(cmdline[9])
    try:
        oid = subprocess.check_output(cmdline).strip()
    except subprocess.CalledProcessError as cpe:
        print("Error compiling {}\n stdout: {}\n stderr: {}".format(
            wf_source,
            cpe.stdout,
            cpe.stderr
        ))
        raise
    return oid.decode("ascii")
    
# TODO run workflow
# def run_workflow(workflow_id):

def test_global_wf_from_wdl():
    # As Alice, compile workflow and create global workflow
    login_alice()

    tname = "global_wf_from_wdl"
    wf_source = os.path.join(os.path.abspath(test_dir), "multi_user", "{}.wdl".format(tname))
    compiler_flags = ["-instanceTypeSelection", "dynamic"]
    workflow_id = build_workflow(
        wf_source=wf_source,
        project_id=PROJECT_ID,
        folder=specific_applet_folder(tname),
        version_id=VERSION_ID,
        compiler_flags=compiler_flags
    )
    full_workflow_id = "{}:{}".format(PROJECT_ID, workflow_id)
    print("Compiled {}".format(full_workflow_id))

    # Strictly increasing global WF version
    global_workflow_version = "1.0.{}".format(int(time.time()))
    print("Global workflow version {}".format(global_workflow_version))

    # Build global workflow from workflow

    # TODO use full_workflow_id after
    # https://jira.internal.dnanexus.com/browse/APPS-975 is fixed
    global_workflow_cmd = [
        "dx",
        "build",
        "--globalworkflow",
        "--publish",
        "--from",
        workflow_id,
        "--version",
        global_workflow_version,
        "--bill-to",
        BILLING_ORG
    ]
    print(" ".join(global_workflow_cmd))

    try:
        subprocess.call(global_workflow_cmd)
    except subprocess.CalledProcessError as cpe:
        print("Error building global workflow from {}\n stdout: {}\n stderr: {}".format(
            workflow_id,
            cpe.stdout,
            cpe.stderr
        ))
        raise

    global_workflow_name = "globalworkflow-{}".format(tname)
    print("Global workflow name {}".format(global_workflow_name))

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

    try:
        subprocess.call(add_developers_cmd)
        subprocess.call(add_users_cmd)
    except subprocess.CalledProcessError as cpe:
        print("Error during developer actions on {}\n stdout: {}\n stderr: {}".format(
            global_workflow_name,
            cpe.stdout,
            cpe.stderr
        ))
        raise

    # TODO test some more developer actions as Alice

    # As Bob, run global workflow
    login_bob()

    # TODO as Bob, run global workflow

    # TODO handle reporting success / failure; use pytest?

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
        "--folder", help="Use an existing folder, instead of building dxCompiler"
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

    version_id = util.get_version_id(top_dir)
    register_version(version_id)

    project = util.get_project(args.project)
    if project is None:
        raise RuntimeError("Could not find project {}".format(args.project))
    print("Project {} ({})".format(project.name, project.get_id()))
    register_project(project.get_id())

    # Do folder setup as Alice
    login_alice()

    if args.folder is None:
        base_folder = util.create_build_dirs(project, version_id)
    else:
        # Use existing prebuilt base folder
        base_folder = args.folder
        util.create_build_subdirs(project, base_folder)
    print("Base folder {}".format(base_folder))
    applet_folder = base_folder + "/applets"
    test_folder = base_folder + "/test"
    register_folders(applet_folder, test_folder)

    # Build the dxCompiler jar file, only on us-east-1
    test_dict = {"aws:us-east-1": project.name + ":" + base_folder}
    assets = util.build(project, base_folder, version_id, top_dir, test_dict,
                        force=False)
    print("assets: {}".format(assets))

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
