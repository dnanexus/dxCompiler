import logging
import os
import shutil
import subprocess as sp
import json
import functools
import time
import dxpy

from glob import glob
from typing import Set, Dict, List
from concurrent import futures
from dxcint.Context import Context
from dxcint.Dependency import Dependency


# enumeration of supported languages
# used in the following ways:
# - language_dir = lang.upper()
# - JAR name = "dxExecutor{}".format(lang) - only first letter is upper case (Cwl, Wdl)
# - asset name = "dx{}rt".format(lang.upper())

def async_retry(max_retries: int = 5, delay: int = 5):
    """
    A decorator method to perform retry a decorated function
    Args:
        max_retries: int. Number of maximum retries.
        delay: int. Amount of time to sleep in seconds between retries.
    Returns: A result of a decorated callable.
    """
    def async_retry_inner(func):
        @functools.wraps(func)
        def async_retry_wrapper(*args, **kwargs):
            lock = args[0].context.lock
            for i in range(0, max_retries):
                try:
                    logging.info(f"Retry {i} for function `{func.__name__}`")
                    ret_value = func(*args, **kwargs)
                    return ret_value
                except Exception:
                    with lock:
                        logging.info(
                            f"Error when running an async retry for function `{func.__name__}`\n"
                            f"With ARGS: {args}\n"
                            f"With KWARGS: {kwargs}\n"
                            f"Retry in {delay} sec"
                        )
                    time.sleep(delay)
            else:
                raise Exception(f"Failed after {max_retries} retries for function `{func.__name__}`\n"
                                f"With ARGS: {args}\n"
                                f"With KWARGS: {kwargs}")

        return async_retry_wrapper

    return async_retry_inner


class TerraformError(Exception):
    """
    Handles Terraform errors
    """


class Terraform(object):
    def __init__(self, languages: Set[str], context: Context, dependencies: Set[Dependency]):
        self._languages = languages
        self._context = context
        self._dependencies = dependencies
        for language in languages:
            self._clean_up(language.upper())
        self._asset_switch = {
            "wdl": self._wdl_asset,
            "cwl": self._cwl_asset,
            "cwl.json": self._cwl_asset
        }
        self._local_created_dirs = {}

    @property
    def context(self):
        return self._context

    def build(self) -> List[str]:
        """
        Main interface for preparing the local space and the platform for testing dxCompiler.

        Returns: List[str]. List of built asset IDs
        """
        build_queue = {self._asset_switch[x] for x in self._languages}     # set of partial functions
        _ = self._generate_config_file()
        _ = self._build_compiler()
        with futures.ThreadPoolExecutor(max_workers=len(build_queue)) as executor:
            future_to_build_task = {executor.submit(build_task): build_task for build_task in build_queue}
            assets = {future_to_build_task[f]: f.result() for f in futures.as_completed(future_to_build_task)}
        return list(assets.values())

    def _clean_up(self, language: str) -> bool:
        """
        Method to remove local asset directories for a given language.
        Args:
            language: str. Language for which to destroy terraformed paths with resources.

        Returns: bool. True if all is done
        """
        language_dir = os.path.join(self._context.repo_root_dir, "applet_resources", language)
        if os.path.exists(language_dir):
            shutil.rmtree(language_dir)
        else:
            logging.info(f"No local asset directories for `{language}` language were pre-existing")
        return True

    def _wdl_asset(self) -> str:
        return self._make_prerequisites("wdl")

    def _cwl_asset(self) -> str:
        return self._make_prerequisites("cwl")

    def _make_prerequisites(self, language: str) -> str:
        always_capital_lang = language.upper()
        language_specific_dependencies = [x for x in self._dependencies if language in x.languages]
        try:
            local_asset_dirs = self._create_local_asset_dir(always_capital_lang)
            for dependency in language_specific_dependencies:
                _ = dependency.link(local_asset_dirs.get("bin"))
                dependency.update_dot_env(local_asset_dirs.get("home"))
            _ = self._create_asset_spec(always_capital_lang)
            asset_id = self._build_asset(always_capital_lang)
            return asset_id
        except Exception as e:
            self._clean_up(always_capital_lang)
            raise e

    @async_retry()
    def _build_asset(self, language: str) -> str:
        asset_name = f"dx{language}rt"
        destination = f"{self._context.project_id}:{self._context.platform_build_dir}/{asset_name}"
        cwd = os.getcwd()
        os.chdir(os.path.join(self._context.repo_root_dir, "applet_resources"))
        logging.info(f"Creating a runtime asset for {language}")
        proc = sp.Popen(
            ["dx", "build_asset", language, "--destination", destination],
            stdout=sp.PIPE,
            stderr=sp.PIPE
        )
        out, err = proc.communicate()
        if proc.returncode != 0:
            raise TerraformError(f"Building DNAnexus asset raised {err.decode()}")
        else:
            with self._context.lock:
                logging.info(f"Building DNAnexus asset returned {out.decode()}\n{err.decode()}")
        os.chdir(cwd)
        asset_id = dxpy.search.find_one_data_object(
            classname="record",
            project=self._context.project_id,
            name=asset_name,
            folder=self._context.platform_build_dir,
            more_ok=False
        )
        return asset_id.get("id")

    def _create_asset_spec(self, language: str) -> Dict:
        spec_exports = [x.export_spec() for x in self._dependencies].remove(None)
        exec_depends = [
                           {"name": "openjdk-8-jre-headless"},
                           {"name": "bzip2"},
                           {"name": "jq"},
                       ] + (spec_exports or [])
        asset_spec = {
            "version": self._context.version,
            "name": f"dx{language}rt",
            "title": f"dx{language} asset",
            "release": "20.04",
            "distribution": "Ubuntu",
            "execDepends": exec_depends,
            "instanceType": "mem1_ssd1_v2_x4",
            "description": f"Prerequisites for running {language} workflows compiled to the platform",
            "excludeResource": ["/dev/console"],
        }
        asset_spec_file = os.path.join(
            self._context.repo_root_dir, "applet_resources", language, "dxasset.json"
        )
        with open(asset_spec_file, "w") as asset_spec_handle:
            asset_spec_handle.write(json.dumps(asset_spec, indent=4))
        return asset_spec

    def _create_local_asset_dir(self, language: str) -> Dict:
        logging.info(f"Creating local asset directories for {language}.")
        language_dir = os.path.join(self._context.repo_root_dir, "applet_resources", language)
        resources_dir = os.path.join(language_dir, "resources")
        local_assets = {
            "resources": resources_dir,
            "bin": os.path.join(resources_dir, "usr", "bin"),
            "home": os.path.join(resources_dir, "home", "dnanexus")
        }
        for key, value in local_assets.items():
            os.makedirs(value, exist_ok=True)
        self._local_created_dirs.update({language: local_assets})
        return local_assets

    def _generate_config_file(self) -> str:
        """
        Create a dxCompiler_runtime.conf file (in typesafe-config format) in the
        compiler's resources directory. It holds a mapping from region to project
        where the runtime asset is stored.

        Returns: str. Config as a string.
        """
        region_project_hocon = []
        all_regions = []
        region = self._context.project_info.get("region")
        project_name = self._context.project_info.get("name")
        dx_path = f"{project_name}:{self._context.platform_build_dir}"
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
        runtime_conf_path = os.path.join(
            self._context.repo_root_dir, "compiler", "src", "main", "resources", "dxCompiler_runtime.conf"
        )
        if os.path.exists(runtime_conf_path):
            os.remove(runtime_conf_path)
        with open(runtime_conf_path, "w") as fd:
            fd.write(conf)
        all_regions_str = ", ".join(all_regions)
        logging.info(
            f"Built configuration regions [{all_regions_str}] into {runtime_conf_path}"
        )
        return conf

    def _build_compiler(self) -> bool:
        try:
            sp.check_call(["sbt", "clean"])
            sp.check_call(["sbt", "assembly"])
        except sp.CalledProcessError as e:
            print(e.stdout)
            print(e.stderr)
            raise e
        jar_exec_origin = os.path.join(
            self._context.repo_root_dir,
            "applet_resources",
            "dxCompiler.jar"
        )
        jar_exec_destination = os.path.join(
            self._context.repo_root_dir,
            f"dxCompiler-{self._context.version}.jar"
        )
        for existing_exe in glob(os.path.join(self._context.repo_root_dir, "dxCompiler-*-SNAPSHOT.jar")):
            os.remove(jar_exec_destination)
            logging.info(f"Removing dxCompiler exe {os.path.basename(existing_exe)}")
        shutil.move(jar_exec_origin, jar_exec_destination)
        return True

