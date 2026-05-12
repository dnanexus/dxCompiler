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
import subprocess
import sys
import os

import dxpy
from dxpy.exceptions import DXJobFailureError

here = os.path.dirname(sys.argv[0])
top_dir = os.path.dirname(os.path.abspath(here))

COPY_FILE_APP_NAME = "dxwdl_copy"
URL_DURATION = 60 * 60 * 24  # 24 hours

HOME_REGION = "aws:us-east-1"
HOME_PROJECT_NAME = "dxCompiler"

# Mapping from region name to destination project name
REGION_TO_PROJECT = {
    "aws:us-east-1":      "dxCompiler",
    "aws:ap-southeast-2": "dxCompiler_Sydney",
    "azure:westus":       "dxCompiler_Azure",
    "azure:westeurope":   "dxCompiler_Amsterdam",
    "aws:eu-central-1":   "dxCompiler_Berlin",
    "aws:eu-west-2":      "dxCompiler_London",
    "aws:eu-west-2-g":    "dxCompiler_Europe_London",
    "azure:uksouth-ofh":  "dxCompiler_OFH_TRE_London",
    "oci:us-ashburn-1":   "dxCompiler_Ashburn",
}


def get_project(project_name):
    """Find a DNAnexus project by name."""
    try:
        project = dxpy.DXProject(project_name)
        return project
    except dxpy.DXError:
        pass
    results = list(dxpy.find_projects(name=project_name, return_handler=True, level="VIEW"))
    if len(results) == 0:
        return None
    if len(results) == 1:
        return results[0]
    # Prefer owned projects
    owned = [r for r in results if r.describe()["level"] == "ADMINISTER"]
    if len(owned) == 1:
        return owned[0]
    raise RuntimeError("Found {} projects named '{}'".format(len(results), project_name))


def find_asset(project, folder, language):
    """Return the AssetBundle record for the given language in the given folder, or None."""
    asset_name = "dx{}rt".format(language.upper())
    assets = list(dxpy.search.find_data_objects(
        classname="record",
        project=project.get_id(),
        name=asset_name,
        folder=folder,
        return_handler=True,
    ))
    if len(assets) == 0:
        return None
    if len(assets) == 1:
        return assets[0]
    raise RuntimeError("Found {} records named '{}' in {}:{}".format(
        len(assets), asset_name, project.name, folder))


def _wait_for_job(job):
    """Block until job completes, printing timestamps every 60 s to keep CI alive."""
    noise = subprocess.Popen(["/bin/bash", "-c", "while true; do sleep 60; date; done"])
    try:
        job.wait_on_done()
    finally:
        noise.kill()


def clone_asset(copy_app, source_record, dest_proj, folder, language):
    """
    Clone the asset from source_record into dest_proj/folder using the dxwdl_copy app.
    Returns the newly created AssetBundle record in the destination project.
    """
    asset_name = "dx{}rt".format(language.upper())

    # Check if the asset already exists in the destination
    existing = find_asset(dest_proj, folder, language)
    if existing is not None:
        print("Asset '{}' already exists in {}:{} ({}). Nothing to do.".format(
            asset_name, dest_proj.name, folder, existing.get_id()))
        return existing

    # Get the underlying file from the source record
    fid = source_record.get_details()['archiveFileId']['$dnanexus_link']
    asset_file_name = dxpy.describe(fid)['name']

    print("Generating pre-authenticated download URL for {} ...".format(asset_file_name))
    url = dxpy.DXFile(fid).get_download_url(
        preauthenticated=True,
        project=dxpy.DXFile.NO_PROJECT_HINT,
        duration=URL_DURATION,
    )[0]

    dest_proj.new_folder(folder, parents=True)
    dest_region = dxpy.describe(dest_proj.get_id())['region']

    print("Launching copy job in region {} (project: {} {}) ...".format(
        dest_region, dest_proj.name, dest_proj.get_id()))
    dxjob = copy_app.run(
        app_input={"url": url, "folder": folder, "filename": asset_file_name},
        name="copy {} to {}".format(asset_name, dest_region),
        project=dest_proj.get_id(),
        priority="high",
    )
    print("Copy job: {}".format(dxjob.get_id()))

    print("Waiting for copy job to complete ...")
    _wait_for_job(dxjob)
    print("Copy job finished.")

    # Find the uploaded file
    results = list(dxpy.find_data_objects(
        classname="file",
        visibility="hidden",
        name=asset_file_name,
        project=dest_proj.get_id(),
        folder=folder,
    ))
    file_ids = [p["id"] for p in results]
    if len(file_ids) == 0:
        raise RuntimeError("Copy job succeeded but no file found at {}:{}/{}".format(
            dest_proj.get_id(), folder, asset_file_name))
    if len(file_ids) > 1:
        raise RuntimeError("Found {} files at {}:{}/{}, expected exactly one".format(
            len(file_ids), dest_proj.get_id(), folder, asset_file_name))

    # Create the AssetBundle record pointing at the copied file
    asset_properties = source_record.get_properties()
    asset_properties['cloned_from'] = source_record.get_id()

    dest_record = dxpy.new_dxrecord(
        name=source_record.name,
        types=['AssetBundle'],
        details={'archiveFileId': dxpy.dxlink(file_ids[0])},
        properties=asset_properties,
        project=dest_proj.get_id(),
        folder=folder,
        close=True,
    )
    print("Created AssetBundle record: {}".format(dest_record.get_id()))
    return dest_record


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

    # Verify the copy app actually supports the destination region
    app_supported_regions = set(copy_app.describe()['regionalOptions'].keys())
    if args.region not in app_supported_regions:
        print("ERROR: Region '{}' is not supported by the '{}' app. "
              "Supported regions: {}".format(
                  args.region, COPY_FILE_APP_NAME,
                  ", ".join(sorted(app_supported_regions))),
              file=sys.stderr)
        sys.exit(1)

    folder = "/releases/{}".format(args.version)
    language = args.language.capitalize()  # "Wdl" or "Cwl"

    # Find the source asset in the home region
    home_proj = get_project(HOME_PROJECT_NAME)
    if home_proj is None:
        raise RuntimeError("Could not find home project '{}'".format(HOME_PROJECT_NAME))
    print("Home project: {} ({})".format(home_proj.name, home_proj.get_id()))

    source_record = find_asset(home_proj, folder, language)
    if source_record is None:
        raise RuntimeError("No {} asset found in {}:{} — has the release been built?".format(
            language, HOME_PROJECT_NAME, folder))
    print("Source asset: {} ({})".format(source_record.name, source_record.get_id()))

    # Find or create destination project
    dest_proj_name = REGION_TO_PROJECT[args.region]
    dest_proj = get_project(dest_proj_name)
    if dest_proj is None:
        raise RuntimeError("Could not find destination project '{}'".format(dest_proj_name))
    print("Destination project: {} ({})".format(dest_proj.name, dest_proj.get_id()))

    clone_asset(copy_app, source_record, dest_proj, folder, language)
    print("Done.")


if __name__ == '__main__':
    main()
