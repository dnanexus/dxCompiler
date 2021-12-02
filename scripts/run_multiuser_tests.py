#!/usr/bin/env python3
import argparse
import dxpy
import os
import sys

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
APPLET_FOLDER = ""
TEST_FOLDER = ""

def register_tokens(alice_token, bob_token):
    global ALICE_SECURITY_CONTEXT
    ALICE_SECURITY_CONTEXT["auth_token"] = alice_token
    global BOB_SECURITY_CONTEXT
    BOB_SECURITY_CONTEXT["auth_token"] = bob_token

def login_alice():
    dxpy.set_api_server_info(host="stagingapi.dnanexus.com")
    dxpy.set_security_context(ALICE_SECURITY_CONTEXT)
    print("Logged in as {}".format(dxpy.bindings.whoami()))

def login_bob():
    dxpy.set_api_server_info(host="stagingapi.dnanexus.com")
    dxpy.set_security_context(BOB_SECURITY_CONTEXT)
    print("Logged in as {}".format(dxpy.bindings.whoami()))

def register_folders(applet_folder, test_folder):
    global APPLET_FOLDER
    APPLET_FOLDER = applet_folder
    global TEST_FOLDER
    TEST_FOLDER = test_folder

def specific_applet_folder(test_name):
    return "{}/{}".format(APPLET_FOLDER, test_name)

# TODO compile workflow
def compile_workflow(test_name, project, version_id, c_flags):
    
# TODO run workflow
# def run_workflow(workflow_id):

def test_global_wf_from_wdl():
    # As Alice, compile workflow and create global workflow
    login_alice()

    # TODO compile workflow

    # TODO determine incremented version

    # TODO make globalworkflow from workflow

    # TODO test developer actions on global workflow

    # As Bob, run global workflow
    login_bob()

    # TODO run global workflow

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

    # Do folder setup as Alice
    login_alice()

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
    
    print("project: {} ({})".format(project.name, project.get_id()))
    print("folder: {}".format(base_folder))
    applet_folder = base_folder + "/applets"
    test_folder = base_folder + "/test"
    register_folders(applet_folder, test_folder)

    test_dict = {"aws:us-east-1": project.name + ":" + base_folder}

    # Build the dxCompiler jar file, only on us-east-1
    assets = util.build(project, base_folder, version_id, top_dir, test_dict,
                        force=args.build is not None)
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
