from dxcint.RegisteredTest import RegisteredTest
import os
import json
from typing import Dict


class ResultsTestMixin(RegisteredTest):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._results_suffix = "_results.json"
        self._results = self._import_results()

    @property
    def results(self) -> Dict:
        return self._results

    def _import_results(self, *args, **kwargs) -> Dict:
        results_basename = f"{self._test_name}{self._results_suffix}"
        results_src = os.path.join(os.path.dirname(self._src_file), results_basename)
        if os.path.exists(results_src):
            with open(results_src, "r") as results_src_handle:
                return json.load(results_src_handle)
        else:
            return {}
