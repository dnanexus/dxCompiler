#!/bin/bash -e

# Create a public docker image for dxCompiler that allows a simple command line invocation

# Find the root directory of the distribution, regardless of where
# the script is called from
base_dir=$(dirname "$0")
root_dir=$base_dir/../..

# Figure out the version number
conf_file=$root_dir/compiler/src/main/resources/application.conf
if [[ ! -e ${conf_file} ]]; then
    echo "Expecting the configuration file to be placed in ${conf_file}"
    exit 1
fi
version=$(grep version ${conf_file} | awk '{print $3}' | sed 's/\"//g' )
echo "dxCompiler version is ${version}"

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

echo "For the next steps to work you need to:"
echo "(1) be logged into docker.io"
echo "(2) have permissions to create a repository for dnanexus"
echo "echo PASSWORD | sudo docker login -u USERNAME --password-stdin"
echo ""
echo "pushing to docker hub"
sudo docker push dnanexus/dxcompiler:${version}
sudo docker push dnanexus/dxcompiler:latest
