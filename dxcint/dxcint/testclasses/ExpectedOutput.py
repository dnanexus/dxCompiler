from typing import Dict, List

from dxcint.RegisteredTest import RegisteredTest
from dxcint.Context import Context
from dxcint.Messenger import State
from dxcint.utils import (
    download_dxfile,
    dict_compare,
    list_dx_folder,
    sort_dicts,
    sort_maybe_mixed,
    async_retry,
    link_to_dxfile,
    get_checksum,
)
from dxcint.mixins.ResultsTestMixin import ResultsTestMixin
import dxpy
import os


class ExpectedOutput(ResultsTestMixin, RegisteredTest):
    def __init__(self, src_file: str, category: str, test_name: str, context: Context):
        super().__init__(src_file, category, test_name, context)
        self._output = None
        self.file_cache = {}
        self.folder_cache = {}

    @async_retry
    def _run_executable(self) -> str:
        execution = self._run_executable_inner()
        return execution.describe().get("id")

    def _extract_outputs(self) -> Dict:
        desc = dxpy.describe(self.job_id)
        if desc["class"] == "analysis":
            if self._locked:
                return desc["output"]
            else:
                stages = desc["stages"]
                for stage in stages:
                    if stage["id"] == "stage-outputs":
                        return stage["execution"]["output"]
                raise RuntimeError(
                    f"Analysis for test {self.name} does not have stage 'outputs'"
                )
        elif desc["class"] == "job":
            return desc["output"]
        else:
            raise RuntimeError(f"Unknown {desc['class']}")

    def _validate_outputs(self, exec_outputs, expected_output) -> bool:
        try:
            if exec_outputs is None:
                if expected_output is None:
                    return True
                else:
                    return False
            for key, expected_val in expected_output:
                exec_name, *field_name_parts = key.split(".")
                field_name1 = ".".join(field_name_parts)
                field_name2 = "___".join(field_name_parts)
                if exec_name != self._test_name:
                    self.context.logger.error(
                        "Execution name in output does not match test name"
                    )
                    return False
                if field_name1 in exec_outputs:
                    exec_val = exec_outputs[field_name1]
                elif field_name2 in exec_outputs:
                    exec_val = exec_outputs[field_name2]
                elif expected_val is None:
                    continue
                else:
                    self.context.logger.error(
                        f"Execution output does not contain expected field: {field_name1}"
                    )
                    return False
                if not self._compare_values(expected_val, exec_val, field_name1):
                    return False

        except Exception:
            return False
        return True

    def _validate(self) -> Dict:
        self.messenger.wait_for_completion()
        if self.messenger.state == State.FINISHED:

            try:
                outputs = self._extract_outputs()
            except RuntimeError as e:
                return {
                    "passed": False,
                    "message": "Extracting outputs threw an exception" + str(e),
                }

            if self._validate_outputs(outputs, self._results):
                return {
                    "passed": True,
                    "message": f"Execution of the test {self.name} passed as expected.",
                }
            else:
                return {
                    "passed": False,
                    "message": f"Results of test {self.name} are invalid.\nExpected: {self._output}.\nReceived: {outputs}",
                }
        else:
            return {
                "passed": False,
                "message": f"Execution of the test {self.name} DID NOT pass as expected.",
            }

    def _unwrap(self, actual, expected):
        if isinstance(actual, dict) and "___" in actual:
            actual = actual["___"]
            if isinstance(expected, dict) and "___" in expected:
                expected = expected["___"]
        if isinstance(actual, dict) and "wrapped___" in actual:
            actual = actual["wrapped___"]
            if isinstance(expected, dict) and "wrapped___" in expected:
                expected = expected["wrapped___"]
        return actual, expected

    def _compare_values(self, expected, actual, field):
        self._unwrap(actual, expected)

        if isinstance(actual, list) and isinstance(expected, list):
            return self._compare_lists(actual, expected, field)
        if isinstance(actual, dict) and isinstance(expected, dict):
            if (
                len(actual) == 1
                and "$dnanexus_link" in actual
                and len(expected) == 1
                and "$dnanexus_link" in expected
            ):
                _, _, modified, _ = dict_compare(actual, expected)
                if modified:
                    self.context.logger.error(
                        f"Given files are not the same ({modified})."
                    )
                return not bool(modified)

        if isinstance(expected, dict) and (
            expected.get("class") in {"File", "Directory"}
            or expected.get("type") in {"File", "Folder"}
        ):
            return self._compare_result_path(actual, expected, field)

        if (
            isinstance(actual, dict)
            and actual.get("type") == "File"
            and "uri" in actual
        ):
            actual = actual["uri"]
        if isinstance(actual, dict) and "$dnanexus_link" in actual:
            actual = download_dxfile(link_to_dxfile(actual, self.context.project_id))

        if isinstance(actual, dict) and isinstance(expected, dict):
            expected_keys = set(expected.keys())
            actual_keys = set(expected.keys())
            if expected_keys != actual_keys:
                self.context.logger.error(
                    f"Analysis {self.name} gave unexpected results"
                )
                self.context.logger.error(
                    f"Field {field} should have keys ({expected_keys}), actual = ({actual_keys})"
                )
                return False
            for k in expected_keys:
                if not self._compare_values(expected[k], actual[k], f"{field}[{k}]"):
                    return False
            else:
                return True

        if str(actual).strip() != str(expected).strip():
            self.context.logger.error(f"Analysis {self.name} gave unexpected results")
            self.context.logger.error(
                f"Field {field} should have keys ({expected_keys}), actual = ({actual_keys})"
            )
            return False

        return True

    def _compare_result_path(self, result, expected_val, field_name: str):
        cls = result.get("class", result.get("type"))
        expected_cls = expected_val.get("class", expected_val.get("type"))
        if cls == "File" and expected_cls == "File":
            return self._compare_result_file(
                result,
                expected_val,
                field_name,
            )
        elif cls in {"Directory", "Folder"} and expected_cls in {"Directory", "Folder"}:
            return self._compare_result_directory(
                result,
                expected_val,
                field_name,
            )
        else:
            self.context.logger.error(f"Analysis {self.name} gave unexpected results")
            self.context.logger.error(
                f"Field {field_name} should be of class ({expected_cls}), actual = ({cls})"
            )
            return False

    def _compare_result_file(self, result, expected_val, field_name: str):
        expected_checksum = (
            expected_val["checksum"] if "checksum" in expected_val else None
        )
        algo = expected_checksum.split("$")[0] if expected_checksum else None

        location = None
        size = None
        checksum = None
        secondary_files = None

        if isinstance(result, str):
            contents = result.strip()
        elif result.get("class", result.get("type")) == "File":
            contents = result.get("contents")
            location = result.get("location", result.get("path", result.get("uri")))
            if isinstance(location, dict) and "$dnanexus_link" in location:
                dxfile = link_to_dxfile(location, self.context.project_id)
                location = os.path.join(
                    dxfile.describe()["folder"], dxfile.describe()["name"]
                )
                if contents is None:
                    contents = download_dxfile(dxfile)
            size = result.get("size")
            checksum = result.get("checksum")
            secondary_files = result.get("secondaryFiles", [])
        elif "$dnanexus_link" in result:
            dxfile = link_to_dxfile(result, self.context.project_id)
            contents = download_dxfile(dxfile)
            location = os.path.join(
                dxfile.describe()["folder"], dxfile.describe()["name"]
            )
        else:
            self.context.logger.error(f"Analysis {self.name} gave unexpected results")
            self.context.logger.error(f"unsupported file value {result}")
            return False

        expected_location = expected_val.get("location", expected_val.get("path"))
        expected_basename = expected_val.get(
            "basename",
            os.path.basename(expected_location) if expected_location else None,
        )
        if expected_basename and expected_basename != "Any":
            basename = os.path.basename(location) if location else None
            if basename != expected_basename:
                self.context.logger.error(
                    f"Analysis {self.name} gave unexpected results"
                )
                self.context.logger.error(
                    f"Field {field_name} should have location/path with basename ({expected_basename}), actual = ({basename})"
                )
                return False

        # the result is a cwl File - match the contents, checksum, and/or size
        if "contents" in expected_val and contents != expected_val["contents"]:
            self.context.logger.error(f"Analysis {self.name} gave unexpected results")
            self.context.logger.error(
                f"Field {field_name} should have contents ({expected_val['contents']}), actual = ({result.get('contents')})"
            )
            return False

        if "size" in expected_val:
            if size is None and contents is not None:
                size = len(contents)
            if size != expected_val["size"]:
                self.context.logger.error(
                    f"Analysis {self.name} gave unexpected results"
                )
                self.context.logger.error(
                    f"Field {field_name} should have size ({expected_val['size']}), actual: ({size})"
                )
                return False

        if expected_checksum:
            checksum = checksum or (get_checksum(contents, algo) if contents else None)
            if checksum != expected_checksum:
                self.context.logger.error(
                    f"Analysis {self.name} gave unexpected results"
                )
                self.context.logger.error(
                    f"Field {field_name} should have checksum ({expected_checksum}), actual = ({checksum})"
                )
                return False

        expected_secondary_files = expected_val.get("secondaryFiles")
        if expected_secondary_files:
            # TODO: sort both lists rather than doing an all-by-all comparison
            for expected in expected_secondary_files:
                for actual in secondary_files:
                    if self._compare_result_path(
                        actual,
                        expected,
                        f"{field_name}.secondaryFiles",
                    ):
                        secondary_files.remove(actual)
                        break
                else:
                    self.context.logger.error(
                        f"Analysis {self.name} gave unexpected results"
                    )
                    self.context.logger.error(
                        f"Field {field_name} is missing secondaryFile ({expected}) from ({secondary_files})"
                    )
                    return False

        return True

    def _compare_lists(self, actual: List, expected: List, field: str) -> bool:
        actual = list(filter(lambda x: x is not None, actual))
        expected = list(filter(lambda x: x is not None, expected))
        n = len(actual)
        if n != len(expected):
            self.context.logger.error(f"Analysis {self.name} gave unexpected results")
            self.context.logger.error(
                f"Field {field} should have length ({len(expected)}), actual = ({len(actual)})"
            )
            return False
        if n == 0:
            return True
        elif n > 1:
            if isinstance(actual[0], dict):
                if (
                    all(
                        len(act) == 1
                        and isinstance(act, dict)
                        and "$dnanexus_link" in act
                        for act in actual
                    )
                ) and (
                    all(
                        len(exp) == 1
                        and isinstance(exp, dict)
                        and "$dnanexus_link" in exp
                        for exp in expected
                    )
                ):
                    for exp, act in zip(expected, actual):
                        if not self._compare_values(exp, act, field):
                            return False
                    return True
                actual, expected = sort_dicts(actual, expected)
            else:
                actual = sort_maybe_mixed(actual)
                expected = sort_maybe_mixed(expected)

        for i, (e, a) in enumerate(zip(expected, actual)):
            if not self._compare_values(e, a, f"{field}[{i}]".format(field, i)):
                return False
        else:
            return True

    def _compare_result_directory(self, result, expected_val, field_name: str) -> bool:
        location = result.get("location", result.get("path", result.get("uri")))
        if location is not None and location.startswith("dx://"):
            project_id, folder = location[5:].split(":")
            project = dxpy.DXProject(project_id)
        else:
            folder = location

        if "basename" in result:
            basename = result["basename"]
        elif folder:
            basename = os.path.basename(folder)
        else:
            basename = None

        expected_location = expected_val.get("location", expected_val.get("path"))
        expected_basename = expected_val.get(
            "basename",
            os.path.basename(expected_location) if expected_location else None,
        )
        if expected_basename and expected_basename != "Any":
            if basename != expected_basename:
                self.context.logger.error(
                    f"Analysis {self.name} gave unexpected results"
                )
                self.context.logger.error(
                    f"Field {field_name} should have location/path with basename ({expected_basename}), actual = ({basename})"
                )
                return False

        expected_listing = expected_val.get("listing")
        if expected_listing:
            if "listing" in result:
                listing = result["listing"]
            elif not folder:
                self.context.logger.error(
                    f"Analysis {self.name} gave unexpected results"
                )
                self.context.logger.error(
                    f"Field {field_name} is missing a folder, actual = ({result})"
                )
                return False
            else:
                listing = list_dx_folder(project, folder)
            listing_len = len(listing) if listing else 0
            if len(expected_listing) != listing_len:
                self.context.logger.error(
                    f"Analysis {self.name} gave unexpected results"
                )
                self.context.logger.error(
                    f"Field {field_name} should have listing ({expected_listing}), actual = ({listing})"
                )
                return False
            for expected in expected_listing:
                for actual in listing:
                    if self._compare_result_path(
                        actual,
                        expected,
                        "{}.listing".format(field_name),
                    ):
                        listing.remove(actual)
                        break
                else:
                    self.context.logger.error(
                        f"Analysis {self.name} gave unexpected results"
                    )
                    self.context.logger.error(
                        f"Field {field_name} is missing item ({expected}) from listing ({listing})"
                    )
                    return False

        return True
