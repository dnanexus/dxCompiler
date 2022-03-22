import logging
import os
import subprocess
from typing import List, Dict
from dxcint.RegisteredTest import RegisteredTest
from dxcint.Context import Context


class Terraform(object):
    def __init__(self, registered_tests: List[RegisteredTest], context: Context):
        self._languages = set(x.language for x in registered_tests)
        self._asset_switch = {
            "wdl": self._wdl_asset,
            "cwl": self._cwl_asset,
            "cwl.json": self._cwl_asset
        }
        self._context = context

    def _build_asset(self, language: str) -> Dict:
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
                logging.info(f"Successfully built asset for language {language}")
        except subprocess.CalledProcessError as e:
            logging.error(e.stdout)
            logging.error(e.stderr)
            raise e
        os.chdir(crnt_work_dir)

    def _wdl_asset(self) -> Dict:
        self._build_asset("wdl")
        pass

    def _cwl_asset(self) -> Dict:
        self._build_asset("cwl")
        pass

    def build(self):
        for language in self._languages:
            self._asset_switch[language]()
