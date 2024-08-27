#!/bin/bash -e

# Specify the dxCompiler version from
# https://github.com/dnanexus/dxCompiler/releases
VERSION=2.11.6

if [[ -z "${VERSION}" ]]; then
  echo "Specify the dxCompiler version in the script"
  exit 1
fi

echo "Downloading dxCompiler v${VERSION} from GitHub"
rm -f dxCompiler-${VERSION}.jar
wget https://github.com/dnanexus/dxCompiler/releases/download/${VERSION}/dxCompiler-${VERSION}.jar

# echo "building a docker image"
# sudo docker build --build-arg VERSION=${version} -t dnanexus/dxcompiler:${version} .

# echo "tagging as latest"
# sudo docker tag dnanexus/dxcompiler:${version} dnanexus/dxcompiler:latest
