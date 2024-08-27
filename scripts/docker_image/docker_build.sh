#!/bin/bash

set -euo pipefail

source "docker.env"

if [[ -z "${VERSION}" ]]; then
  echo "Specify the dxCompiler version in docker.env"
  exit 1
fi

echo "Downloading dxCompiler v${VERSION} from GitHub"
rm -f dxCompiler-${VERSION}.jar
wget https://github.com/dnanexus/dxCompiler/releases/download/${VERSION}/dxCompiler-${VERSION}.jar

echo "Building Docker image"
docker build --build-arg VERSION=${VERSION} -t dnanexus/dxcompiler:${VERSION} .
docker tag dnanexus/dxcompiler:${VERSION} dnanexus/dxcompiler:latest
