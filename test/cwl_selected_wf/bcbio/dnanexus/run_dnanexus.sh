#!/bin/bash
set -eu -o pipefail

#DX_PROJECT_ID=brad_cwl_gvcf
DX_PROJECT_ID=project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq # dxCompiler-playground
PNAME=gvcf-joint

dx mkdir -p /$PNAME-workflow
dx upload $PNAME-workflow/main-$PNAME-samples.json --path $DX_PROJECT_ID:/$PNAME-workflow/
bcbiovm_python ~/drive/az/nexus/dx-cwl/dx-cwl run-workflow /dx-cwl-run/main-$PNAME/main-$PNAME /$PNAME-workflow/main-$PNAME-samples.json --project $DX_PROJECT_ID --token $DX_AUTH_TOKEN
