#!/bin/bash
set -eo pipefail

# List of test project ids
TEST_PROJECTS=('project-G2054Z80139VJGzy0Kq02PZ1')

# Check user is logged in
username=$(dx whoami --id)
rc=$?
if [[ $rc -ne 0 ]]; then
  echo "Could not get username; you are probably not logged in to DNAnexus"
  exit $rc
fi

cd test-dxcompiler-in-project
pwd

# Build, run test applet in each project
for p in ${TEST_PROJECTS[@]}; do
  dx select $p
  echo "Building and running test applet in $p"
  dx build -f --run --priority normal -y
done
