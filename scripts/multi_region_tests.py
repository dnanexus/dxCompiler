#!/usr/bin/env python3
from __future__ import print_function

import argparse
import dxpy
import json
import pprint
import os
import time
import sys
import subprocess
import util
from dxpy.exceptions import DXJobFailureError

######################################################################
# multi-region test.
#
# Compile the trivial workflow on all supported regions, and see that it runs.

here = os.path.dirname(sys.argv[0])
top_dir = os.path.dirname(os.path.abspath(here))
test_dir = os.path.join(os.path.abspath(top_dir), "test")

# Load region-to-project mapping from the shared config file.
# To add or retire a region, edit scripts/regions.json instead of this file.
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "regions.json")) as _f:
    REGION_TO_PROJECT = json.load(_f)["regions"]

# Ordered list of project names derived from the mapping
projects = list(REGION_TO_PROJECT.values())

target_folder = "/release_test"

def wait_for_completion(test_exec_objs):
    print("awaiting completion ...")
    # wait for analysis to finish while working around Travis 10m console inactivity timeout
    noise = subprocess.Popen(["/bin/bash", "-c", "while true; do sleep 60; date; done"])
    try:
        for exec_obj in test_exec_objs:
            exec_obj.wait_on_done()
    finally:
        noise.kill()
    print("done")

# Run [workflow] on several inputs, return the analysis ID.
def run_workflow(dx_proj, test_folder, oid):
    dx_proj.new_folder(test_folder, parents=True)
    wf = dxpy.DXWorkflow(project=dx_proj.get_id(), dxid=oid.decode("utf-8"))
    return wf.run({},
                  project=dx_proj.get_id(),
                  folder=test_folder)

# Build a workflow.
#
# wf             workflow name
# classpath      java classpath needed for running compilation
# folder         destination folder on the platform
def build_test(source_file, dx_proj, folder, version_id):
    dx_proj.new_folder(folder, parents=True)
    print("Compiling {} to project {}:/{}".format(source_file, dx_proj.name, folder))
    cmdline = [ "java", "-jar",
                os.path.join(top_dir, "dxCompiler-{}.jar".format(version_id)),
                "compile",
                source_file,
                "-force",
                "-locked",
                "-folder", folder,
                "-project", dx_proj.get_id() ]
    print(" ".join(cmdline))
    oid = subprocess.check_output(cmdline).strip()
    return oid

def main():
    argparser = argparse.ArgumentParser(description="Run compiler tests on the platform")
    argparser.add_argument("--compile-only", help="Only compile the workflows, don't run them",
                           action="store_true", default=False)
    argparser.add_argument("--skip-regions",
                           help="Comma-separated regions to skip testing (e.g. aws:eu-west-2). "
                                "Uses the same region names as --skip-clone-regions in build_release.py.",
                           default="")
    args = argparser.parse_args()

    # strip whitespace around each entry to tolerate "aws:eu-west-2, azure:westeurope"
    skip_regions = (
        set(r.strip() for r in args.skip_regions.split(",") if r.strip())
        if args.skip_regions else set()
    )

    # Validate region names against known mapping
    if skip_regions:
        unknown = skip_regions - set(REGION_TO_PROJECT.keys())
        if unknown:
            print("ERROR: Unknown region(s) in --skip-regions: {}. "
                  "Known regions: {}".format(", ".join(sorted(unknown)),
                                             ", ".join(sorted(REGION_TO_PROJECT.keys()))))
            sys.exit(1)

    # Build a reverse map: project_name -> region, to filter out skipped regions
    project_to_region = {v: k for k, v in REGION_TO_PROJECT.items()}
    active_projects = [p for p in projects if project_to_region.get(p) not in skip_regions]

    if skip_regions:
        skipped = [p for p in projects if p not in active_projects]
        print("Skipping projects for regions {}: {}".format(skip_regions, skipped))

    version_id = util.get_version_id(top_dir)
    wdl_source_file = os.path.join(test_dir, "multi_region/trivial.wdl")
    dx_objects = []
    test_exec_objs=[]

    # build version of the applet on all regions
    for proj_name in active_projects:
        dx_proj = util.get_project(proj_name)
        if dx_proj is None:
            raise RuntimeError("Could not find project {}".format(proj_name))
        oid = build_test(wdl_source_file, dx_proj, target_folder, version_id)
        if args.compile_only:
            continue
        anl = run_workflow(dx_proj, target_folder, oid)
        print("Running {}".format(oid))
        test_exec_objs.append(anl)

    # Wait for completion
    print("executables: " + ", ".join([a.get_id() for a in test_exec_objs]))
    wait_for_completion(test_exec_objs)

if __name__ == '__main__':
    main()
