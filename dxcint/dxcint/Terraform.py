import logging
import os
import shutil
import subprocess
import json
import functools
import time
import dxpy

from typing import Set, Dict, Optional, List
from concurrent import futures
from dxcint.Context import Context
from dxcint.Dependency import Dependency


# enumeration of supported languages
# used in the following ways:
# - language_dir = lang.upper()
# - JAR name = "dxExecutor{}".format(lang)
# - asset name = "dx{}rt".format(lang.upper())

class Terraform(object):
    def __init__(self, languages: Optional[Set[str]], context: Context, dependencies: List[Dependency]):
        self._languages = languages
        self._asset_switch = {
            "wdl": self._wdl_asset,
            "cwl": self._cwl_asset,
            "cwl.json": self._cwl_asset
        }
        self._context = context
        self._dependencies = dependencies
        self._local_created_dirs = {}

    def build(self) -> List[str]:
        """
        Main interface for preparing the local space and the platform for testing dxCompiler.
        :return: List[str]. List of built asset IDs
        """
        build_queue = set([self._asset_switch[x] for x in self._languages])     # set of partial functions
        self._generate_config_file()
        with futures.ThreadPoolExecutor(max_workers=len(build_queue)) as executor:
            future_to_build_task = {executor.submit(build_task): build_task for build_task in build_queue}
            assets = {future_to_build_task[f]: f.result() for f in futures.as_completed(future_to_build_task)}
        return list(assets.values())

    def destroy(self, language: str) -> bool:
        """
        Method to remove local asset directories for a given language.
        :param language:
        :return: bool. True if all is done
        """
        local_asset_dirs = self._local_created_dirs.pop(language, None)
        if local_asset_dirs:
            for key, value in local_asset_dirs.items():
                if os.path.exists(value):
                    logging.info(f"Removing local asset directory {value}")
                    shutil.rmtree(value)
                else:
                    continue
        else:
            logging.info(f"No local asset directories were created for `{language.upper()}` language")
        return True

    def _wdl_asset(self) -> str:
        return self._make_prerequisites("wdl")

    def _cwl_asset(self) -> str:
        return self._make_prerequisites("cwl")

    def _async_retry(self, max_retries: int = 5, delay: int = 5):
        """
        A decorator method to perform retry a decorated function
        :param max_retries: int.
        :param delay: int. Amount of time to sleep in seconds between retries.
        :return: A result of a decorated callable.
        """
        def _async_retry_inner(func):
            @functools.wraps(func)
            def _async_retry_wrapper(*args, **kwargs):
                for i in range(0, max_retries):
                    try:
                        logging.info(f"Retry {i} for function `{func.__name__}`")
                        func(*args, **kwargs)
                        break
                    except Exception:
                        with self._context.lock:
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
            return _async_retry_wrapper
        return _async_retry_inner

    def _make_prerequisites(self, language: str):
        language_specific_dependencies = [x for x in self._dependencies if language in x.languages]
        try:
            local_asset_dirs = self._create_local_asset_dir(language)
            for dependency in language_specific_dependencies:
                _ = dependency.link(local_asset_dirs.get("bin"))
                dependency.update_dot_env(local_asset_dirs.get("home"))
            _ = self._create_asset_spec(language)
            asset_id = self._build_asset(language)
            return asset_id
        except Exception as e:
            self.destroy(language)
            raise e

    @_async_retry()
    def _build_asset(self, language: str) -> str:
        asset_name = f"dx{language.upper()}rt"
        destination = f"{self._context.project_id}:{self._context.platform_build_dir}/{asset_name}"
        cwd = os.getcwd()
        os.chdir(os.path.join(self._context.repo_root_dir, "applet_resources"))
        logging.info(f"Creating a runtime asset for {language}")
        try:
            subprocess.run(
                ["dx", "build_asset", language.upper(), "--destination", destination],
                check=True,
                universal_newlines=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            with self._context.lock:
                logging.info(f"Successfully built asset for language {language}")
        except subprocess.CalledProcessError as e:
            logging.error(e.stdout)
            logging.error(e.stderr)
            raise e
        os.chdir(cwd)
        return dxpy.search.find_one_data_object(
            classname="record",
            project=self._context.project_id,
            name=asset_name,
            folder=self._context.platform_build_dir,
            more_ok=False
        )

    def _create_asset_spec(self, language: str) -> Dict:
        exec_depends = [
                           {"name": "openjdk-8-jre-headless"},
                           {"name": "bzip2"},
                           {"name": "jq"},
                       ] + (self._dependencies or [])
        asset_spec = {
            "version": self._context.version,
            "name": f"dx{language.upper()}rt",
            "title": f"dx{language.upper()} asset",
            "release": "20.04",
            "distribution": "Ubuntu",
            "execDepends": exec_depends,
            "instanceType": "mem1_ssd1_v2_x4",
            "description": f"Prerequisites for running {language.upper()} workflows compiled to the platform",
            "excludeResource": ["/dev/console"],
        }
        asset_spec_file = os.path.join(
            self._context.repo_root_dir, "applet_resources", language.upper(), "dxasset.json"
        )
        with open(asset_spec_file, "w") as asset_spec_handle:
            asset_spec_handle.write(json.dumps(asset_spec, indent=4))
        return asset_spec

    def _create_local_asset_dir(self, language: str) -> Dict:
        logging.info(f"Creating local asset directories for {language}.")
        language_dir = os.path.join(self._context.repo_root_dir, "applet_resources", language.upper())
        if os.path.exists(language_dir):
            shutil.rmtree(language_dir)
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

    def _generate_config_file(self) -> None:
        """
        Create a dxCompiler_runtime.conf file (in typesafe-config format) in the
        compiler's resources directory. It holds a mapping from region to project
        where the runtime asset is stored.
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

        rt_conf_path = os.path.join(
            self._context.repo_root_dir, "compiler", "src", "main", "resources", "dxCompiler_runtime.conf"
        )
        if os.path.exists(rt_conf_path):
            os.remove(rt_conf_path)
        with open(rt_conf_path, "w") as fd:
            fd.write(conf)
        all_regions_str = ", ".join(all_regions)
        logging.info(
            f"Built configuration regions [{all_regions_str}] into {rt_conf_path}"
        )
