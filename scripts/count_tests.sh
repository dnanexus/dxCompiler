#!/bin/bash
set -e -u -o pipefail
cd $1
tools_count=0
run_count=0
for f in $(ls -1 *.cwl.json); do
  ((tools_count=tools_count+1))
  rc=$((ls -la "${f/.cwl.json/_input}"[0-9]".json" 2>/dev/null || echo "1")|wc -l)
  ((run_count=run_count+rc))
done
echo "$1 tool_test_count $tools_count run_count $run_count"