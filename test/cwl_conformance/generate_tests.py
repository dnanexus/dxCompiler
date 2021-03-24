# this script processes inputs.yml and generates the input.json and results.json
# files used by the integration test script
# TODO: right now this script expects the cwl-v1.2 repo to be checked out at the
#  same level as the dxCompiler project, and for this script to be run from within
#  the test/cwl_conformance/tools dir. It should be re-written to read all the files
#  directly from the repo https://github.com/common-workflow-language/cwl-v1.2
import json
import os
from pathlib import Path
import shutil
import yaml

project_id = "project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq"


def convert_files(raw_value):
    if isinstance(raw_value, dict):
        if "class" in raw_value:
            cls = raw_value["class"]
            if cls == "File":
                if "location" in raw_value:
                    raw_value["location"] = "dx://{}:/test_data/cwl/{}".format(
                        project_id, raw_value["location"]
                    )
                return raw_value
            else:
                return dict(
                    (k, convert_files(v))
                    for k, v in raw_value.items()
                )
        else:
            return raw_value
    elif isinstance(raw_value, list):
        return [convert_files(x) for x in raw_value]
    else:
        return raw_value

def create_test(tool_file, test, index, tools_dir):
    filename = os.path.splitext(tool_file.name)[0]
    basedir = test["basedir"]

    if test.get("job"):
        job_path = basedir / test["job"]
        with open(job_path) as inp:
            if job_path.suffix in (".yml", ".yaml"):
                raw_inputs = yaml.load(inp)
            else:
                raw_inputs = json.load(inp)
        inputs = dict(
            ("{}.{}".format(filename, arg_name), convert_files(arg_value))
            for (arg_name, arg_value) in raw_inputs.items()
        )
        if index < 0:
            input_name = "{}_input.json".format(filename)
        else:
            input_name = "{}_input{}.json".format(filename, index)
        with open(tools_dir / input_name, "w") as out:
            json.dump(inputs, out, indent=2)

    if "output" in test:
        if index < 0:
            output_name = "{}_results.json".format(filename)
        else:
            output_name = "{}_results{}.json".format(filename, index)
        results = dict(
            ("{}.{}".format(filename, arg_name), arg_value)
            for arg_name, arg_value in test["output"].items()
        )
        with open(tools_dir / output_name, "w") as out:
            json.dump(results, out, indent=2)


def generate_tests(tests_path, cwl_path, tools_dir):
    with open(tests_path) as inp:
        tests = yaml.load(inp)
    test_dict = {}
    for test in tests:
        if "$import" in test:
            generate_tests(cwl_path / test["$import"], cwl_path, tools_dir)
        else:
            tool_path = tests_path.parent / test["tool"]
            if not os.path.exists(tools_dir / tool_path.name):
                continue
            if tool_path in test_dict:
                test_list = test_dict[tool_path]
            else:
                test_list = []
                test_dict[tool_path] = test_list
            test["basedir"] = tests_path.parent
            test_list.append(test)

    for tool_file, test_list in test_dict.items():
        if len(test_list) == 1:
            create_test(tool_file, test_list[0], -1, tools_dir)
        else:
            for i, test in enumerate(test_list, 1):
                create_test(tool_file, test, i, tools_dir)


if __name__ == "__main__":
    cwl_path = Path("../../../cwl-v1.2")
    tools_dir = Path("workflows")
    generate_tests(cwl_path / "conformance_tests.yaml", cwl_path, tools_dir)
