#!/usr/bin/env python3
# This script builds a dxCompiler release. A release consists of the client (dxCompiler.jar),
# which is uploaded to

import argparse
import dxpy
import json
import os
import subprocess
import time
import sys

import util

here = os.path.dirname(sys.argv[0])
top_dir = os.path.dirname(os.path.abspath(here))

HOME_REGION = "aws:us-east-1"
URL_DURATION = 60 * 60 * 24
SLEEP_TIME = 5
COPY_FILE_APP_NAME = "dxwdl_copy"
COPY_FILE_APP = dxpy.find_one_app(zero_ok=False, more_ok=False, name=COPY_FILE_APP_NAME, return_handler=True)
NUM_RETRIES = 1

TEST_DICT = {
    "aws:us-east-1" :  "dxCompiler_playground"
}

# Load region-to-project mapping from the shared config file.
# To add or retire a region, edit scripts/regions.json instead of this file.
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "regions.json")) as _f:
    RELEASE_DICT = json.load(_f)["regions"]



def main():
    argparser = argparse.ArgumentParser(description="Build a dxCompiler release")
    argparser.add_argument("--force",
                           help="Build even if there is an existing version",
                           action='store_true',
                           default=False)
    argparser.add_argument("--multi-region",
                           help="Copy to all supported regions",
                           action='store_true',
                           default=False)
    argparser.add_argument("--dry-run",
                           help="Don't build any artifacts",
                           action='store_true',
                           default=False)
    argparser.add_argument("--skip-clone-regions",
                           help="Comma-separated regions to skip asset cloning (e.g. aws:eu-west-2). "
                                "Region is still included in JAR config.",
                           default="")
    args = argparser.parse_args()

    # build multi-region jar for releases, or
    # if explicitly specified
    multi_region = args.multi_region
    # strip whitespace around each entry to tolerate "aws:eu-west-2, azure:westeurope"
    skip_clone_regions = (
        set(r.strip() for r in args.skip_clone_regions.split(",") if r.strip())
        if args.skip_clone_regions else set()
    )

    # Choose which dictionary to use
    if multi_region:
        project_dict = RELEASE_DICT
    else:
        project_dict = TEST_DICT

    # Validate skip_clone_regions against known regions and disallow home region
    if skip_clone_regions:
        unknown = skip_clone_regions - set(project_dict.keys())
        if unknown:
            print("ERROR: Unknown region(s) in --skip-clone-regions: {}. "
                  "Known regions: {}".format(", ".join(sorted(unknown)),
                                             ", ".join(sorted(project_dict.keys()))),
                  file=sys.stderr)
            sys.exit(1)
        if HOME_REGION in skip_clone_regions:
            print("ERROR: Cannot skip the home region ({}).".format(HOME_REGION),
                  file=sys.stderr)
            sys.exit(1)

    project = util.get_project(project_dict[HOME_REGION])
    print("project: {} ({})".format(project.name, project.get_id()))

    # Figure out what the current version is
    version_id = util.get_version_id(top_dir)
    print("version: {}".format(version_id))

    # Set the folder
    folder = "/releases/{}".format(version_id)
    print("folder: {}".format(folder))

    if args.dry_run:
        args.force = False

    # remove the existing directory paths
    if args.force:
        for proj_name in project_dict.values():
            print("removing path {}:{}".format(proj_name, folder))
            dx_proj = util.get_project(proj_name)
            try:
                dx_proj.remove_folder(folder, recurse=True)
            except dxpy.DXError:
                pass

    # Make sure the target directory exists
    project.new_folder(folder, parents=True)

    # Build the asset, and the compiler jar file.
    path_dict = dict(map(lambda kv: (kv[0], kv[1] + ":" + folder),
                         project_dict.items()))
    if args.dry_run:
        return

    home_ad = util.build(project, folder, version_id, top_dir, path_dict)

    if multi_region:
        for lang, asset_desc in home_ad.items():
            home_rec = dxpy.DXRecord(asset_desc.asset_id)
            all_regions = project_dict.keys()

            # Leave only regions where the asset is missing and not explicitly skipped
            target_regions = []
            for dest_region in all_regions:
                if dest_region in skip_clone_regions:
                    print("Skipping asset clone for region: {}".format(dest_region), file=sys.stderr)
                    continue
                dest_proj = util.get_project(project_dict[dest_region])
                dest_asset = util.find_asset(dest_proj, folder, lang)
                if dest_asset == None:
                    target_regions.append(dest_region)

            util.clone_asset(COPY_FILE_APP, home_rec, folder, target_regions, project_dict)

if __name__ == '__main__':
    main()
