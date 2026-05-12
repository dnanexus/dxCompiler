#!/usr/bin/env python3
"""
Standalone script to copy a dxCompiler runtime asset from the home region
(aws:us-east-1) into a specific destination region.

Use this to backfill a region that was skipped during a release
(e.g. because it was unavailable at release time).

Usage:
    python3 copy_asset_to_region.py \
        --version 2.x.y \
        --region aws:eu-west-2 \
        --language wdl \
        [--token <dx-auth-token>]

If --token is not supplied the script uses the currently active dx login session.
"""
from __future__ import print_function

import argparse
import json
import sys
import os

import dxpy
import util

COPY_FILE_APP_NAME = "dxwdl_copy"

HOME_REGION = "aws:us-east-1"
HOME_PROJECT_NAME = "dxCompiler"

# Load region-to-project mapping from the shared config file.
# To add or retire a region, edit scripts/regions.json instead of this file.
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "regions.json")) as _f:
    REGION_TO_PROJECT = json.load(_f)["regions"]


def main():
    argparser = argparse.ArgumentParser(
        description="Copy a dxCompiler runtime asset from the home region to a destination region"
    )
    argparser.add_argument("--version", required=True,
                           help="dxCompiler release version (e.g. 2.11.0)")
    argparser.add_argument("--region", required=True,
                           help="Destination region (e.g. aws:eu-west-2)")
    argparser.add_argument("--language", required=True, choices=["wdl", "cwl"],
                           help="Runtime language asset to copy")
    argparser.add_argument("--token",
                           help="DNAnexus auth token. If omitted, uses the current dx session.")
    args = argparser.parse_args()

    # Authenticate
    if args.token:
        dxpy.set_security_context({"auth_token_type": "Bearer", "auth_token": args.token})

    # Validate region
    if args.region not in REGION_TO_PROJECT:
        print("ERROR: Unknown region '{}'. Supported regions: {}".format(
            args.region, ", ".join(sorted(REGION_TO_PROJECT.keys()))))
        sys.exit(1)
    if args.region == HOME_REGION:
        print("ERROR: Destination region is the home region ({}). Nothing to copy.".format(HOME_REGION))
        sys.exit(1)

    # Initialize copy app AFTER login (avoids module-level import-time failure)
    print("Looking up copy app '{}' ...".format(COPY_FILE_APP_NAME))
    copy_app = dxpy.find_one_app(
        zero_ok=False, more_ok=False, name=COPY_FILE_APP_NAME, return_handler=True
    )
    print("Found copy app: {}".format(copy_app.get_id()))

    folder = "/releases/{}".format(args.version)
    language = args.language.capitalize()  # "Wdl" or "Cwl"

    # Find the source asset in the home region
    home_proj = util.get_project(HOME_PROJECT_NAME)
    if home_proj is None:
        raise RuntimeError("Could not find home project '{}'".format(HOME_PROJECT_NAME))
    print("Home project: {} ({})".format(home_proj.name, home_proj.get_id()))

    source_record = util.find_asset(home_proj, folder, language)
    if source_record is None:
        raise RuntimeError("No {} asset found in {}:{} — has the release been built?".format(
            language, HOME_PROJECT_NAME, folder))
    print("Source asset: {} ({})".format(source_record.name, source_record.get_id()))

    # Find destination project
    dest_proj_name = REGION_TO_PROJECT[args.region]
    dest_proj = util.get_project(dest_proj_name)
    if dest_proj is None:
        raise RuntimeError("Could not find destination project '{}'".format(dest_proj_name))
    print("Destination project: {} ({})".format(dest_proj.name, dest_proj.get_id()))

    # Idempotency check: skip if AssetBundle record already exists
    existing = util.find_asset(dest_proj, folder, language)
    if existing is not None:
        print("Asset 'dx{}rt' already exists in {}:{} ({}). Nothing to do.".format(
            language.upper(), dest_proj.name, folder, existing.get_id()))
        return

    util.clone_asset(copy_app, source_record, folder, [args.region], REGION_TO_PROJECT)
    print("Done.")


if __name__ == '__main__':
    main()
