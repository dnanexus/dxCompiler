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

echo "Downloading dxCompiler v${VERSION} from GitHub"
rm -f "${base_dir}/dxCompiler-${VERSION}.jar"
wget https://github.com/dnanexus/dxCompiler/releases/download/${VERSION}/dxCompiler-${VERSION}.jar -O "${base_dir}/dxCompiler-${VERSION}.jar"

# echo "Building Docker image"
docker build -f "${base_dir}/Dockerfile" --build-arg VERSION=${VERSION} --build-arg BASE_DIR=${base_dir} -t dnanexus/dxcompiler:${VERSION} .
