#!/bin/bash -e

# TODO Specify the version

# TODO Download dxCompiler-${version}.jar

if [[ ! -f dxCompiler-${version}.jar ]]; then
    # To build the docker image, we need a copy of the jar file. We
    # download it from github to make sure we have an up to date version.
    rm -f dxCompiler-${version}.jar
    wget https://github.com/dnanexus/dxCompiler/releases/download/${version}/dxCompiler-${version}.jar
fi

echo "building a docker image"
sudo docker build --build-arg VERSION=${version} -t dnanexus/dxcompiler:${version} .

echo "tagging as latest"
sudo docker tag dnanexus/dxcompiler:${version} dnanexus/dxcompiler:latest
