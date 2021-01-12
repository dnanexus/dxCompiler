#!/bin/bash

# Requires bash v4

declare -A PROJECTS
PROJECTS[dxCompiler]=aws:us-east-1
PROJECTS[dxCompiler_Berlin]=aws:eu-central-1
PROJECTS[dxCompiler_Amsterdam]=azure:westeurope
PROJECTS[dxCompiler_Sydney]=aws:ap-southeast-2
PROJECTS[dxCompiler_London]=aws:eu-west-2

for i in "${!PROJECTS[@]}"; do
    proj_id=$(dx new project $i --region ${PROJECTS[$i]} --bill-to=org-dnanexus_apps --brief)
    echo $i $proj_id
done
