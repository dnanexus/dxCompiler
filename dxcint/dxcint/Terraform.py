import logging
import os
import shutil
import subprocess
import json
import sys
import time
from typing import Set, Dict, Optional, List
from dxcint.Context import Context
from dxcint.Dependency import Dependency


# enumeration of supported languages
# used in the following ways:
# - language_dir = lang.upper()
# - subproject = "executor{}".format(lang)
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

        self._local_asset_dirs = self._create_local_asset_dir(language)

    def destroy(self) -> bool:
        """
        Method to remove local asset directories
        :return: bool. True, when all is done
        """
        for key, value in self._local_asset_dirs.items():
            if os.path.exists(value):
                logging.info(f"Removing local asset directory {value}")
                shutil.rmtree(value)
            else:
                continue
        return True

    def _make_prerequisites(
            self,
            language,
            lock
    ):

        # Link in the shared resources
        for dependency in self._dependencies:
            dependency.link(self._local_asset_dirs.get("bin"))

        # Link in executor-specific resources, if any
        # TODO
        source_lang_res_dir = os.path.join(
            self._context.repo_root_dir, f"executor{language}", "applet_resources"
        )
        if os.path.exists(source_lang_res_dir):
            for f in os.listdir(source_lang_res_dir):
                os.link(os.path.join(source_lang_res_dir, f), os.path.join(language_dir, f))

        # Create the asset description file
        self._create_asset_spec(language)

        # Create the .env file if necessary
        if env_vars:
            self._create_env_file()

        # Create an asset from the executor jar file and its dependencies,
        # this speeds up applet creation.
        destination = f"{self._context.project_id}:{self._context.platform_build_dir}/dx{language.upper()}rt"
        for i in range(0, max_num_retries):
            try:
                with lock:
                    logging.info(f"Creating a runtime asset for {language} (try {i})")
                self._build_asset(top_dir, language, destination, lock)
                break
            except Exception:
                with lock:
                    logging.info(
                        f"Error creating runtime asset for {language}; sleeping for 5 seconds before trying again",
                        sys.exc_info()
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

    def _build_asset(self, language: str) -> Dict:
        cwd = os.getcwd()
        os.chdir(os.path.join(self._context.repo_root_dir, "applet_resources"))
        try:
            subprocess.run(
                ["dx", "build_asset", language.upper(), "--destination", destination],
                check=True,
                universal_newlines=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            with lock:
                logging.info(f"Successfully built asset for language {language}")
        except subprocess.CalledProcessError as e:
            logging.error(e.stdout)
            logging.error(e.stderr)
            raise e
        os.chdir(cwd)

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
        asset_spec_file = os.path.join(self._context.repo_root_dir, "applet_resources", language.upper(), "dxasset.json")
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
        return local_assets

    def _create_env_file(self):
        dot_env = "\n".join(f"{key}={val}" for key, val in env_vars.items())
        dot_env_file = os.path.join(self._local_asset_dirs.get("home"), ".env")
        with open(dot_env_file, "wt") as dot_env_handle:
            dot_env_handle.write(dot_env)

    def _wdl_asset(self) -> Dict:
        self._build_asset("wdl")
        pass

    def _cwl_asset(self) -> Dict:
        self._build_asset("cwl")
        pass

    def build(self):
        try:
            for language in self._languages:
                self._create_local_asset_dir(language)
                self._asset_switch[language]()
        except Exception as e:
            self.destroy()
            raise e
