import json
import os

from dx_cwl_runner import utils
from dx_cwl_runner.dx import Dx


def get_modified_output(output: dict, outdir, dx: Dx):
    file_id = next(iter(output.values()))
    output_path = dx.get_file_name(outdir, file_id)
    utils.run_cmd(f"dx download {next(iter(output.values()))} --output '{output_path}' --no-progress")
    return {
        "class": "File",
        "checksum": str(utils.get_checksum(output_path)),
        "location": output_path,
        "size": os.path.getsize(output_path)
    }


def create_dx_output(job, process_file, outdir, dx: Dx):
    description = json.loads(utils.run_cmd(f"dx describe {job} --json"))
    outputs = description.get("output")
    results = {}
    if description.get("state") == "done":
        if outputs:
            for output in outputs:
                if type(outputs[output]) is list:
                    results[f"{process_file}.{output}"] = []
                    for idx, value in enumerate(outputs[output]):
                        results[f"{process_file}.{output}"][idx] = get_modified_output(
                            outputs[outputs][idx], outdir, dx=dx
                        )
                if type(outputs[output]) is dict:
                    results[f"{process_file}.{output}"] = get_modified_output(
                        outputs[output], outdir, dx=dx
                    )
                else:
                    results[f"{process_file}.{output}"] = outputs[output]
                    continue
    return results


def write_output(file_name, results, print_output=True):
    with open(file_name, 'w') as result_file:
        json.dump(results, result_file)
    if print_output:
        print(json.dumps(results, indent=4))
