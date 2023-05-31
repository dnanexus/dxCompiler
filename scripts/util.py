#!/usr/bin/env python3
from collections import namedtuple
from concurrent import futures
import dxpy
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
import threading
import traceback

AssetDesc = namedtuple("AssetDesc", "region asset_id project")

max_num_retries = 5
# enumeration of supported languages
# used in the following ways:
# - language_dir = lang.upper()
# - subproject = "executor{}".format(lang)
# - JAR name = "dxExecutor{}".format(lang)
# - asset name = "dx{}rt".format(lang.upper())
languages = ["Wdl", "Cwl"]

DEFAULT_INSTANCE_TYPE = "mem1_ssd1_v2_x4"
GIT_REVISION = subprocess.check_output(
    ["git", "describe", "--always", "--dirty", "--tags"]
).strip()


def info(msg, ex=None):
    print(msg, file=sys.stderr)
    if ex:
        traceback.print_exception(*ex)


# Extract version_id from configuration file
def get_version_id(top_dir):
    appl_conf_path = os.path.join(
        top_dir, "core", "src", "main", "resources", "application.conf"
    )
    pattern = re.compile(r"^(\s*)(version)(\s*)(=)(\s*)(\S+)(\s*)$")
    with open(appl_conf_path, "r") as fd:
        for line in fd:
            line_clean = line.replace('"', "").replace("'", "")
            m = re.match(pattern, line_clean)
            if m is not None:
                return m.group(6).strip()
    raise Exception("version ID not found in {}".format(appl_conf_path))


def find_asset(project, folder, language):
    # get asset_id
    asset_name = "dx{}rt".format(language.upper())
    assets = list(
        dxpy.search.find_data_objects(
            classname="record",
            project=project.get_id(),
            name=asset_name,
            folder=folder,
            return_handler=True,
        )
    )
    if len(assets) == 0:
        return None
    if len(assets) == 1:
        return assets[0]
    raise Exception(
        "More than one asset with name {} found in {}:{}".format(
            asset_name, project, folder
        )
    )


def create_build_subdirs(project, base_folder):
    """Creates subfolder in the base folder needed for running tests"""
    applet_folder = base_folder + "/applets"
    test_folder = base_folder + "/test"
    project.new_folder(test_folder, parents=True)
    project.new_folder(applet_folder, parents=True)


def create_build_dirs(project, version_id):
    username = dxpy.whoami()
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

    project = list(
        dxpy.find_projects(name=project_name, return_handler=True, level="VIEW")
    )
    if len(project) == 0:
        info("Did not find project {0}".format(project_name))
        return None
    elif len(project) == 1:
        return project[0]
    else:
        raise Exception("Found more than 1 project matching {0}".format(project_name))


# download dxda for linux, and place it in the resources
# sub-directory.
def _download_dxda_into_resources(top_dir, dxda_version):
    # TODO: if dxda_version is None, fetch latest
    dxda_dir = os.path.join(top_dir, "applet_resources", "dxda", dxda_version)
    dxda_exe = os.path.join(dxda_dir, "dx-download-agent")

    if not os.path.exists(dxda_exe):
        os.makedirs(dxda_dir, exist_ok=True)
        os.chdir(dxda_dir)
        try:
            # download dxda release, and place it in the resources directory
            if dxda_version.startswith("v"):
                # A proper download-agent release, it starts with a "v"

                subprocess.check_call(
                    [
                        "wget",
                        "https://github.com/dnanexus/dxda/releases/download/{}/dx-download-agent-linux".format(
                            dxda_version
                        ),
                        "-O",
                        "dx-download-agent",
                    ]
                )
            else:
                # A snapshot of the download-agent development branch
                snapshot_dxda_tar = "resources/dx-download-agent-linux.tar"
                command = """sudo docker run --rm --entrypoint=\'\' dnanexus/dxda:{} cat /builds/dx-download-agent-linux.tar > {}""".format(
                    dxda_version, snapshot_dxda_tar
                )
                p = subprocess.Popen(
                    command,
                    universal_newlines=True,
                    shell=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                )
                text = p.stdout.read()
                retcode = p.wait()
                info("downloading dxda {} {}".format(retcode, text))
                subprocess.check_call(
                    ["tar", "-C", "resources", "-xvf", snapshot_dxda_tar]
                )
                os.rename(
                    "resources/dx-download-agent-linux/dx-download-agent",
                    "dx-download-agent",
                )
                os.remove(snapshot_dxda_tar)
                shutil.rmtree("resources/dx-download-agent-linux")
        except subprocess.CalledProcessError as e:
            print(e.stdout)
            print(e.stderr)
            raise e
        # make sure the binary is executable
        os.chmod("dx-download-agent", 0o775)

    return dxda_exe


def _download_dxfuse_to_resources(top_dir, dxfuse_version):
    # TODO: if dxfuse_version is None, fetch latest
    dxfuse_dir = os.path.join(top_dir, "applet_resources", "dxfuse", dxfuse_version)
    dxfuse_exe = os.path.join(dxfuse_dir, "dxfuse")

    if not os.path.exists(dxfuse_exe):
        os.makedirs(dxfuse_dir, exist_ok=True)
        info("downloading dxfuse {} to {}".format(dxfuse_version, dxfuse_exe))
        try:
            subprocess.check_call(
                [
                    "wget",
                    "https://github.com/dnanexus/dxfuse/releases/download/{}/dxfuse-linux".format(
                        dxfuse_version
                    ),
                    "-O",
                    dxfuse_exe,
                ]
            )
        except subprocess.CalledProcessError as e:
            print(e.stdout)
            print(e.stderr)
            raise e

        os.chmod(dxfuse_exe, 0o775)

    return dxfuse_exe


def _download_awscli_to_resources(top_dir, awscli_version):
    awscli_dir = os.path.join(top_dir, "applet_resources", "awscli", awscli_version)
    awscli_zip = os.path.join(awscli_dir, "awscliv2.zip")

    if not (os.path.exists(awscli_zip)):
        os.makedirs(awscli_dir, exist_ok=True)
        info("downloading awscli {} to {}".format(awscli_version, awscli_zip))
        try:
            subprocess.check_call(
                [
                    "wget",
                    "https://awscli.amazonaws.com/awscli-exe-linux-x86_64-{}.zip".format(
                        awscli_version
                    ),
                    "-O",
                    awscli_zip,
                ]
            )
        except subprocess.CalledProcessError as e:
            print(e.stdout)
            print(e.stderr)
            raise e

    return awscli_zip


def _create_asset_spec(version_id, top_dir, language, dependencies=None):
    exec_depends = [
        {"name": "openjdk-11-jre-headless"},
        {"name": "bzip2"},
        {"name": "jq"},
    ] + (dependencies or [])
    asset_spec = {
        "version": version_id,
        "name": "dx{}rt".format(language.upper()),
        "title": "dx{} asset".format(language.upper()),
        "release": "20.04",
        "distribution": "Ubuntu",
        "execDepends": exec_depends,
        "instanceType": "mem1_ssd1_v2_x4",
        "description": "Prerequisites for running {} workflows compiled to the platform".format(
            language.upper()
        ),
        "excludeResource": ["/dev/console"],
    }
    with open(
        os.path.join(top_dir, "applet_resources", language.upper(), "dxasset.json"), "w"
    ) as fd:
        fd.write(json.dumps(asset_spec, indent=4))


# Build a dx-asset from the runtime library.
# Go to the applet_resources directory, before running "dx"
def _build_asset(top_dir, language, destination, lock):
    crnt_work_dir = os.getcwd()
    # build the platform asset
    os.chdir(os.path.join(os.path.abspath(top_dir), "applet_resources"))
    try:
        subprocess.run(
            ["dx", "build_asset", language.upper(), "--destination", destination],
            check=True,
            universal_newlines=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        with lock:
            print("Successfully built asset for language {}".format(language))
            # print(proc.stdout)
            # print(proc.stderr)
    except subprocess.CalledProcessError as e:
        print(e.stdout)
        print(e.stderr)
        raise e
    os.chdir(crnt_work_dir)


def _make_prerequisites(
    project,
    folder,
    version_id,
    top_dir,
    language,
    resources,
    lock,
    dependencies=None,
    env_vars=None,
):
    # Create a folder for the language-specific asset
    language_dir = os.path.join(top_dir, "applet_resources", language.upper())
    build_lang_res_dir = os.path.join(language_dir, "resources")
    build_lang_bin_dir = os.path.join(build_lang_res_dir, "usr", "bin")
    os.makedirs(build_lang_bin_dir, exist_ok=True)

    # Link in the shared resources
    for res in resources:
        os.link(res, os.path.join(build_lang_bin_dir, Path(res).name))

    # Link in executor-specific resources, if any
    source_lang_res_dir = os.path.join(
        top_dir, "executor{}".format(language), "applet_resources"
    )
    if os.path.exists(source_lang_res_dir):
        for f in os.listdir(source_lang_res_dir):
            os.link(os.path.join(source_lang_res_dir, f), os.path.join(language_dir, f))

    # Create the asset description file
    _create_asset_spec(version_id, top_dir, language, dependencies)

    # Create the .env file if necessary
    if env_vars:
        dot_env = "\n".join("{}={}".format(key, val) for key, val in env_vars.items())
        # files in home dir are not included in the final asset
        build_lang_home_dir = os.path.join(build_lang_res_dir, "home", "dnanexus")
        os.makedirs(build_lang_home_dir, exist_ok=True)
        dot_env_file = os.path.join(build_lang_home_dir, ".env")
        with open(dot_env_file, "wt") as out:
            out.write(dot_env)

    # Create an asset from the executor jar file and its dependencies,
    # this speeds up applet creation.
    destination = "{}:{}/dx{}rt".format(project.get_id(), folder, language.upper())
    for i in range(0, max_num_retries):
        try:
            with lock:
                info("Creating a runtime asset for {} (try {})".format(language, i))
            _build_asset(top_dir, language, destination, lock)
            break
        except Exception:
            with lock:
                info(
                    f"Error creating runtime asset for {language}; sleeping for 5 seconds before trying again",
                    sys.exc_info(),
                )
            time.sleep(5)
    else:
        raise Exception("Failed to build the {} runtime asset".format(language))

    # make sure the asset exists and is findable
    asset = find_asset(project, folder, language)
    if asset is None:
        raise Exception(
            "unable to discover the asset created at {}".format(destination)
        )
    return asset


# Create a dxCompiler_runtime.conf file (in typesafe-config format) in the
# compiler's resources directory. It holds a mapping from region to project
# where the runtime asset is stored.
def _gen_config_file(top_dir, project_dict):
    # Create a record for each region
    region_project_hocon = []
    all_regions = []
    for region, dx_path in project_dict.items():
        record = "\n".join(
            [
                "    {",
                '      region = "{}"'.format(region),
                '      path = "{}"'.format(dx_path),
                "    }",
            ]
        )
        region_project_hocon.append(record)
        all_regions.append(region)

    buf = "\n".join(region_project_hocon)
    conf = "\n".join(
        ["dxCompiler {", "  regionToProject = [\n{}\n  ]".format(buf), "}"]
    )

    rt_conf_path = os.path.join(
        top_dir, "compiler", "src", "main", "resources", "dxCompiler_runtime.conf"
    )
    if os.path.exists(rt_conf_path):
        os.remove(rt_conf_path)
    with open(rt_conf_path, "w") as fd:
        fd.write(conf)
    all_regions_str = ", ".join(all_regions)
    info(
        "Built configuration regions [{}] into {}".format(all_regions_str, rt_conf_path)
    )


# Build a fat jar file using sbt-assembly.
# Make sure to run with the working directory being the top dir of the project.
def _sbt_assembly(top_dir):
    os.chdir(os.path.abspath(top_dir))
    jar_paths = dict(
        (
            "dxExecutor{}".format(lang),
            os.path.join(
                top_dir,
                "applet_resources",
                lang.upper(),
                "resources",
                "dxExecutor{}.jar".format(lang),
            ),
        )
        for lang in languages
    )
    jar_paths.update(
        {"dxCompiler": os.path.join(top_dir, "applet_resources", "dxCompiler.jar")}
    )
    for jar_path in jar_paths.values():
        if os.path.exists(jar_path):
            os.remove(jar_path)
    try:
        subprocess.check_call(["sbt", "clean"])
        subprocess.check_call(["sbt", "assembly"])
    except subprocess.CalledProcessError as e:
        print(e.stdout)
        print(e.stderr)
        raise e
    for jar_path in jar_paths.values():
        if not os.path.exists(jar_path):
            raise Exception(f"sbt assembly failed: {jar_path} not found")
    return jar_paths


def build(
    project, folder, version_id, top_dir, path_dict, dependencies=None, force=False
):
    assets = dict((lang, find_asset(project, folder, lang)) for lang in languages)
    build_assets = force or not all(assets.values())

    if build_assets:
        for lang in languages:
            language_dir = os.path.join(top_dir, "applet_resources", lang.upper())
            if os.path.exists(language_dir):
                shutil.rmtree(language_dir)

        if dependencies is None:
            with open(
                os.path.join(top_dir, "scripts/bundled_dependencies.json"), "rt"
            ) as inp:
                dependencies = json.load(inp)

        # get a copy of the download agent (dxda)
        dxda_exe = _download_dxda_into_resources(top_dir, dependencies["dxda"])
        # get a copy of the dxfuse executable
        dxfuse_exe = _download_dxfuse_to_resources(top_dir, dependencies["dxfuse"])
        # get a copy of awscliv2
        resources = [dxda_exe, dxfuse_exe]

        # Create a configuration file
        _gen_config_file(top_dir, path_dict)
        jar_paths = _sbt_assembly(top_dir)
        info("jar_paths: {}".format(jar_paths))

        exec_depends = dependencies.get("execDepends", {})
        env_vars = dependencies.get("env", {})
        # build assets in parallel
        lock = threading.Lock()
        with futures.ThreadPoolExecutor(max_workers=len(languages)) as executor:
            future_to_lang = dict(
                (
                    executor.submit(
                        _make_prerequisites,
                        project,
                        folder,
                        version_id,
                        top_dir,
                        lang,
                        resources,
                        lock,
                        exec_depends.get(lang.lower()),
                        env_vars.get(lang.lower()),
                    ),
                    lang,
                )
                for lang in languages
            )
            assets = dict(
                (future_to_lang[f], f.result())
                for f in futures.as_completed(future_to_lang)
            )

        for (prefix, jar_path) in jar_paths.items():
            # Move the file to the top level directory
            all_in_one_jar = os.path.join(
                top_dir, "{}-{}.jar".format(prefix, version_id)
            )
            shutil.move(jar_path, all_in_one_jar)

        # delete the language-specific dirs
        for lang in languages:
            shutil.rmtree(os.path.join(top_dir, "applet_resources", lang.upper()))

    region = dxpy.describe(project.get_id())["region"]
    asset_descs = dict(
        (lang, AssetDesc(region, asset.get_id(), project))
        for (lang, asset) in assets.items()
    )

    return asset_descs


# Read a JSON file
def read_json_file(path):
    with open(path, "r") as fd:
        return json.load(fd)


# Verify a JSON file
def verify_json_file(path):
    try:
        read_json_file(path)
    except Exception:
        raise RuntimeError("Error verifying JSON file {}".format(path))


# Build an executable
#
# return            Id of the compiled executable
#
# source_file       Workflow source file
# project           Destination project on platform
# folder            Destination folder on platform
# top_dir           Local folder containing dxCompiler.jar
# version_id        dxCompiler version
# compiler_flags    Additional dxCompiler flags
def build_executable(
    source_file, project, folder, top_dir, version_id, compiler_flags=None
):
    cmdline = [
        "java",
        "-jar",
        os.path.join(top_dir, "dxCompiler-{}.jar".format(version_id)),
        "compile",
        source_file,
        "-force",
        "-folder",
        folder,
        "-project",
        project.get_id(),
    ]
    cmdline += compiler_flags or []
    try:
        print(" ".join(cmdline))
        oid = subprocess.check_output(cmdline).strip()
    except subprocess.CalledProcessError as cpe:
        print(
            "Error compiling {}\nstdout: {}\nstderr: {}".format(
                source_file, cpe.stdout, cpe.stderr
            )
        )
        raise
    return oid.decode("ascii")


# Run an executable on 0-many inputs.
#
# return                        List of tuples (input #, analysis or job id)
#
# oid                           Id of executable to run, or name for global workflow
# project                       Destination project on platform
# test_folder                   Destination folder on platform
# test_name                     Test name
# test_inputs                   Inputs for running, if non-default
# debug_flag                    Keep jobs open for debugging?
# delay_workspace_destruction   Delay workspace destruction?
# instance_type                 Instance type, if non-default
# expected_failures             List of test names where failure is expected
def run_executable(
    oid,
    project,
    test_folder,
    test_name,
    test_inputs=None,
    debug_flag=False,
    delay_workspace_destruction=False,
    instance_type=DEFAULT_INSTANCE_TYPE,
    expected_failures=None,
):
    test_inputs = test_inputs or []
    expected_failures = expected_failures or set()

    def once(i):
        try:
            if len(test_inputs) == 0 or i < 0:
                print("Running with empty input")
                inputs = {}
            else:
                print("Running with input file: {}".format(test_inputs[i]))
                inputs = read_json_file(test_inputs[i])
            project.new_folder(test_folder, parents=True)
            if "globalworkflow-" in oid:
                global_workflow_name = oid.replace("globalworkflow-", "")
                exec_obj = dxpy.DXGlobalWorkflow(name=global_workflow_name)
                run_kwargs = {"ignore_reuse_stages": ["*"]}
            elif "workflow-" in oid:
                exec_obj = dxpy.DXWorkflow(project=project.get_id(), dxid=oid)
                run_kwargs = {"ignore_reuse_stages": ["*"]}
            elif "applet-" in oid:
                exec_obj = dxpy.DXApplet(project=project.get_id(), dxid=oid)
                run_kwargs = {"ignore_reuse": True}
            else:
                raise RuntimeError(
                    "Unable to determine executable kind for {}".format(oid)
                )

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
                name="{} {}".format(test_name, GIT_REVISION),
                **run_kwargs,
            )
        except Exception as e:
            print("Exception message {}".format(e))
            return None

    def run(i):
        for _ in range(1, 5):
            retval = once(i)
            if retval is not None:
                return retval
            print("Sleeping for 5 seconds before trying again")
            time.sleep(5)
        else:
            if (
                test_name in expected_failures
                or "{}.{}".format(test_name, i) in expected_failures
            ):
                print("Analysis {}.{} failed as expected".format(test_name, i))
                return None
            else:
                raise RuntimeError("Error running workflow {}".format(oid))

    n = len(test_inputs)
    if n == 0:
        return [(0, run(-1))]
    else:
        return [(i, run(i)) for i in range(n)]
