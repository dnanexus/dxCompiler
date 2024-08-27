#!/bin/bash

set -euo pipefail

# Specify the dxCompiler version from
# https://github.com/dnanexus/dxCompiler/releases
VERSION=2.11.6

if [[ -z "${VERSION}" ]]; then
  echo "Specify VERSION = the dxCompiler version to use"
  exit 1
fi

eval "$(dx env --bash)"

docker run \
     --rm \
     -e DX_SECURITY_CONTEXT="$DX_SECURITY_CONTEXT" \
     -e DX_APISERVER_PROTOCOL="$DX_APISERVER_PROTOCOL" \
     -e DX_APISERVER_HOST="$DX_APISERVER_HOST" \
     -e DX_APISERVER_PORT="$DX_APISERVER_PORT" \
     -e DX_PROJECT_CONTEXT_ID="$DX_PROJECT_CONTEXT_ID" \
     --volume "$PWD":/workdir \
     --workdir /workdir \
     dnanexus/dxcompiler:"$VERSION" "$@"
