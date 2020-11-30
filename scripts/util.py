#!/usr/bin/env python
from __future__ import print_function
from collections import namedtuple
import dxpy
import json
import pwd
import os
import re
import shutil
import subprocess
import sys
import time

AssetDesc = namedtuple('AssetDesc', 'region asset_id project')

max_num_retries = 5
# enumeration of supported languages
# used in the following ways:
# - language_dir = lang.upper()
# - subproject = "executor{}".format(lang)
# - JAR name = "dxExecutor{}".format(lang)
# - asset name = "dx{}rt".format(lang.upper())
languages = ["Wdl"]


# Extract version_id from configuration file
def get_version_id(top_dir):
    appl_conf_path = os.path.join(top_dir, "core", "src", "main", "resources", "application.conf")
    pattern = re.compile(r"^(\s*)(version)(\s*)(=)(\s*)(\S+)(\s*)$")
    with open(appl_conf_path, 'r') as fd:
        for line in fd:
            line_clean = line.replace("\"", "").replace("'", "")
            m = re.match(pattern, line_clean)
            if m is not None:
                return m.group(6).strip()
    raise Exception("version ID not found in {}".format(appl_conf_path))


def find_asset(project, folder, language):
    # get asset_id
    asset_name = "dx{}rt".format(language.upper())
    assets = list(dxpy.search.find_data_objects(classname="record",
                                                project=project.get_id(),
                                                name=asset_name,
                                                folder=folder,
                                                return_handler=True))
    if len(assets) == 0:
        return None
    if len(assets) == 1:
        return assets[0]
    raise Exception("More than one asset with name {} found in {}:{}".format(asset_name, project, folder))


def create_build_subdirs(project, base_folder):
    """Creates subfolder in the base folder needed for running tests"""
    applet_folder = base_folder + "/applets"
    test_folder = base_folder + "/test"
    project.new_folder(test_folder, parents=True)
    project.new_folder(applet_folder, parents=True)


def create_build_dirs(project, version_id):
    user_desc = pwd.getpwuid(os.getuid())
    username = user_desc.pw_name
    base_folder = "/builds/{}/{}".format(username, version_id)
    create_build_subdirs(project, base_folder)
    return base_folder


def get_project(project_name):
    """Try to find the project with the given name or id."""

    # First, see if the project is a project-id.
    try:
        project = dxpy.DXProject(project_name)
        return project
    except dxpy.DXError:
        pass

    project = dxpy.find_projects(name=project_name, return_handler=True, level="VIEW")
    project = [p for p in project]
    if len(project) == 0:
        print('Did not find project {0}'.format(project_name), file=sys.stderr)
        return None
    elif len(project) == 1:
        return project[0]
    else:
        raise Exception('Found more than 1 project matching {0}'.format(project_name))


# download dxda for linux, and place it in the resources
# sub-directory.
def _download_dxda_into_resources(top_dir, dxda_version):
    # TODO: if dxda is already downloaded, check that it's version matches
    # TODO: if dxda_version is None, fetch latest
    os.chdir(os.path.join(top_dir, "applet_resources"))

    # download dxda release, and place it in the resources directory
    if dxda_version.startswith("v"):
        # A proper download-agent release, it starts with a "v"
        subprocess.check_call([
            "wget",
            "https://github.com/dnanexus/dxda/releases/download/{}/dx-download-agent-linux".format(dxda_version),
            "-O",
            "resources/usr/bin/dx-download-agent"])
    else:
        trg_dxda_tar = "resources/dx-download-agent-linux.tar"
        # A snapshot of the download-agent development branch
        command = """sudo  docker run --rm --entrypoint=\'\' dnanexus/dxda:{} cat /builds/dx-download-agent-linux.tar > {}""".format(
            dxda_version, trg_dxda_tar)
        p = subprocess.Popen(command, universal_newlines=True, shell=True,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        text = p.stdout.read()
        retcode = p.wait()
        print("downloading dxda{} {}".format(retcode, text))
        subprocess.check_call(["tar", "-C", "resources", "-xvf", trg_dxda_tar])
        os.rename("resources/dx-download-agent-linux/dx-download-agent",
                  "resources/usr/bin/dx-download-agent")
        os.remove(trg_dxda_tar)
        shutil.rmtree("resources/dx-download-agent-linux")

    os.chmod("resources/usr/bin/dx-download-agent", 0o775)


def _download_dxfuse_to_resources(top_dir, dxfuse_version):
    # TODO: if dxfuse is already downloaded, check that it's version matches
    # TODO: if dxfuse_version is None, fetch latest
    p = os.path.join(top_dir, "applet_resources/resources/usr/bin/dxfuse")
    if not os.path.exists(p):
        print("downloading dxfuse {} to {}".format(dxfuse_version, p))
        subprocess.check_call([
            "wget",
            "https://github.com/dnanexus/dxfuse/releases/download/{}/dxfuse-linux".format(dxfuse_version),
            "-O",
            p])
    os.chmod(p, 0o775)


def _create_asset_spec(version_id, top_dir, language):
    asset_spec = {
        "version": version_id,
        "name": "dx{}rt".format(language.upper()),
        "title": "dx {} asset".format(language.upper()),
        "release": "20.04",
        "distribution": "Ubuntu",
        "execDepends": [
            {"name": "openjdk-8-jre-headless"},
            {"name": "bzip2"},
            {"name": "jq"}
        ],
        "instanceType": "mem1_ssd1_v2_x4",
        "description": "Prerequisites for running {} workflows compiled to the platform".format(language.upper())
    }
    with open(os.path.join(top_dir, "applet_resources", language, "dxasset.json"), 'w') as fd:
        fd.write(json.dumps(asset_spec, indent=4))


# Build a dx-asset from the runtime library.
# Go to the applet_resources directory, before running "dx"
def _build_asset(top_dir, language, destination):
    crnt_work_dir = os.getcwd()
    # build the platform asset
    os.chdir(os.path.join(os.path.abspath(top_dir), "applet_resources"))
    subprocess.check_call(["dx", "build_asset", language, "--destination", destination])
    os.chdir(crnt_work_dir)


def _make_prerequisites(project, folder, version_id, top_dir, language):
    # Create a folder for the language-specific asset
    language_resources_dir = os.path.join(top_dir, "applet_resources", language, "resources")
    os.makedirs(language_resources_dir)

    # Link in the shared resources
    resources = os.path.join(top_dir, "applet_resources", "resources", "*")
    subprocess.check_call(["ln", "-s", resources, language_resources_dir])

    # Create the asset description file
    _create_asset_spec(version_id, top_dir, language)

    # Create an asset from the dxWDL jar file and its dependencies,
    # this speeds up applet creation.
    destination = "{}:{}/dx{}rt".format(project.get_id(), folder, language.upper())
    for i in range(0, max_num_retries):
        try:
            print("Creating a runtime asset for {} (try {})".format(language, i))
            _build_asset(top_dir, language, destination)
            break
        except:
            print("Sleeping for 5 seconds before trying again")
            time.sleep(5)
    else:
        raise Exception("Failed to build the {} runtime asset".format(language))


# Create a dxCompiler_runtime.conf file (in typesafe-config format) in the
# compiler's resources directory. It holds a mapping from region to project
# where the runtime asset is stored.
def _gen_config_file(top_dir, project_dict):
    # Create a record for each region
    region_project_hocon = []
    all_regions = []
    for region, dx_path in project_dict.items():
        record = "\n".join(["    {",
                            '      region = "{}"'.format(region),
                            '      path = "{}"'.format(dx_path),
                            "    }"])
        region_project_hocon.append(record)
        all_regions.append(region)

    buf = "\n".join(region_project_hocon)
    conf = "\n".join(["dxCompiler {",
                      "  regionToproject = [\n{}\n  ]".format(buf),
                      "}"])

    rt_conf_path = os.path.join(top_dir, "compiler", "src", "main", "resources", "dxCompiler_runtime.conf")
    if os.path.exists(rt_conf_path):
        os.remove(rt_conf_path)
    with open(rt_conf_path, 'w') as fd:
        fd.write(conf)
    all_regions_str = ", ".join(all_regions)
    print("Built configuration regions [{}] into {}".format(all_regions_str, rt_conf_path))


# Build a fat jar file using sbt-assembly.
# Make sure to run with the working directory being the top dir of the project.
def _sbt_assembly(top_dir):
    os.chdir(os.path.abspath(top_dir))
    jar_paths = dict(
        ("dxExecutor{}".format(lang), os.path.join(
            top_dir, "applet_resources", lang.upper(), "resources", "dxExecutor{}.jar".format(lang)
        )) for lang in languages.items()
    )
    jar_paths.update({"dxCompiler": os.path.join(top_dir, "applet_resources", "dxCompiler.jar")})
    for jar_path in jar_paths.values():
        if os.path.exists(jar_path):
            os.remove(jar_path)
    subprocess.check_call(["sbt", "clean"])
    subprocess.check_call(["sbt", "assembly"])
    for jar_path in jar_paths.values():
        if not os.path.exists(jar_path):
            raise Exception("sbt assembly failed")
    return jar_paths


def build(project, folder, version_id, top_dir, path_dict, dependencies=None, force=False):
    build_assets = force or not all(find_asset(project, folder, lang) for lang in languages)

    if build_assets:
        if dependencies is None:
            with open(os.path.join(top_dir, "scripts/bundled_dependencies.json"), "rt") as inp:
                dependencies = json.load(inp)

        # make sure the resources directory exists
        if not os.path.exists(os.path.join(top_dir, "applet_resources/resources/usr/bin")):
            os.makedirs(os.path.join(top_dir, "applet_resources/resources/usr/bin"))

        # get a copy of the dxfuse executable
        _download_dxfuse_to_resources(top_dir, dependencies["dxfuse"])

        # get a copy of the download agent (dxda)
        _download_dxda_into_resources(top_dir, dependencies["dxda"])

        _make_prerequisites(project, folder, version_id, top_dir)
        asset = find_asset(project, folder)

        # Create a configuration file
        _gen_config_file(top_dir, path_dict)
        jar_paths = _sbt_assembly(top_dir)

        for (prefix, jar_path) in jar_paths.items():
            # Move the file to the top level directory
            all_in_one_jar = os.path.join(top_dir, "{}-{}.jar".format(prefix, version_id))
            shutil.move(jar_path, all_in_one_jar)

    region = dxpy.describe(project.get_id())['region']
    ad = AssetDesc(region, asset.get_id(), project)

    # Hygiene, remove the new configuration file, we
    # don't want it to leak into the next build cycle.
    # os.remove(crnt_conf_path)
    return ad
