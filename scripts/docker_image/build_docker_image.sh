#!/bin/bash -e

# Specify the dxCompiler version from
# https://github.com/dnanexus/dxCompiler/releases
VERSION=2.11.6

if [[ -z "${VERSION}" ]]; then
  echo "Specify VERSION = the dxCompiler version to use"
  exit 1
fi

echo "Downloading dxCompiler v${VERSION} from GitHub"
rm -f dxCompiler-${VERSION}.jar
wget https://github.com/dnanexus/dxCompiler/releases/download/${VERSION}/dxCompiler-${VERSION}.jar

echo "Building Docker image"
docker build --build-arg VERSION=${VERSION} -t dnanexus/dxcompiler:${VERSION} .
docker tag dnanexus/dxcompiler:${VERSION} dnanexus/dxcompiler:latest
