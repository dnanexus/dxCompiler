#!/bin/bash
set -eu -o pipefail

# dx env
#DX_PROJECT_ID=brad_cwl_gvcf
DX_PROJECT_ID=project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq # dxCompiler-playground
PNAME=gvcf-joint

bcbiovm_python ~/bio/dx-cwl/dx-cwl compile-workflow $PNAME-workflow/main-$PNAME.cwl --project $DX_PROJECT_ID --token $DX_AUTH_TOKEN --assets record-F8YXfF80BV9Pb2G7P4x9fvB5
