#!/bin/bash

set -euo pipefail

# Find the root directory of the distribution, regardless of where
# the script is called from
base_dir=$(dirname "$0")
source "${base_dir}/docker.env"

if [[ -z "${VERSION}" ]]; then
  echo "Specify the dxCompiler version in docker.env"
  exit 1
fi

# bash "${base_dir}/docker_build.sh" 2>&1

output=$(bash "${base_dir}/docker_run.sh" version 2>&1)
if ! echo "$output" | grep -q "${VERSION}"
then
  echo "Expected dxCompiler version ${VERSION}, output was ${output}"
  exit 1
fi
