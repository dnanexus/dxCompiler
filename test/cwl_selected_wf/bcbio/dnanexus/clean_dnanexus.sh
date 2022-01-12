#!/bin/bash
set -eu -o pipefail

rm -rf dx-cwl-run
dx rm -r .cwl_workflow_archive/
dx rm -r .cwl_workflow_output_archive/
dx rm -r main-gvcf-joint-output/
dx rm -r gvcf-joint-workflow/
dx rm -r main-gvcf-joint/
dx rm -r wf-alignment/
dx rm -r wf-variantcall/
dx rm -r wf-jointcall/
dx rm -r dx-cwl-run/
dx rm -r md5sum/
dx rm -r md5sum-output/
dx rm -r wf-alignment-output/
dx rm -r wf-variantcall-output/
dx rm -r wf-jointcall-output/
dx ls
