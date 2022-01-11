#!/bin/bash
set -eu -o pipefail

PNAME=gvcf-joint

rm -rf $PNAME-workflow
bcbio_vm.py cwl --systemconfig=bcbio_system-dnanexus.yaml $PNAME.yaml

