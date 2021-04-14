import dxpy
import json
import utils
import os
from datetime import datetime
def get_env_variables():
    with open(f"{os.environ.get('HOME')}/.dnanexus_config/environment.json") as config:
        env = json.load(config)
        current_dx_folder = env.get("DX_CLI_WD", "/")
    current_dx_project = dxpy.PROJECT_CONTEXT_ID
    input_dx_folder = f"{current_dx_folder}/cwl_runner_{datetime.now().strftime('%Y.%d.%m_%H-%M-%S')}"
    utils.run_cmd(f"dx mkdir -p {current_dx_project}:{input_dx_folder}")
    return current_dx_folder, current_dx_project, input_dx_folder

current_dx_folder, current_dx_project, test_dx_folder = get_env_variables()
