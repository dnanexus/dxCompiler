#!/bin/bash

set -euo pipefail

# This tests that the Dockerfile produces a valid Docker image with Java and
# dxCompiler installed, based on the dxCompiler version specified in docker.env
# (not necessarily the latest version being released, as it isn't yet published).

# Find the root directory of the distribution, regardless of where
# the script is called from
base_dir=$(dirname "$0")
source "${base_dir}/docker.env"

if [[ -z "${VERSION}" ]]; then
  echo "Specify the dxCompiler version in docker.env"
  exit 1
fi

if ! bash "${base_dir}/docker_build.sh" 2>&1
then
  echo "Docker image build failed"
  exit 1
fi

output=$(bash "${base_dir}/docker_run.sh" version 2>&1)
echo "Output from running dxCompiler Docker image to check version:"
echo "${output}"
if ! echo "${output}" | grep -q "${VERSION}"
then
  echo "Expected dxCompiler version ${VERSION}, output was ${output}"
  exit 1
fi
