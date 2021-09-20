#!/bin/bash

# Requirements: base64, jq, gzip

# Usage:
# ./extract_source_code.sh workflow-xxxx
# ./extract_source_code.sh applet-yyyy

set -exo pipefail

exec_id=$1

dx describe $exec_id --json \
| jq .details.sourceCode \
| awk '{$1=$1;print}' \
| sed 's/\"//g' \
| base64 --decode \
| gunzip
