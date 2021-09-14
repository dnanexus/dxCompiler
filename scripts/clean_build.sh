#!/bin/bash

TEST=${1}

# Clean artifacts from prev. builds, local
sbt clean && sbt cleanFiles
find . -name target | xargs rm -rf
rm -rf applet_resources
rm dx*.jar

# Clean artifacts from prev. builds, platform
username=$(dx whoami --id)
rc=$?
if [[ $rc -ne 0 ]]; then
  echo "Could not get username; you are probably not logged in to DNAnexus"
  exit $rc
fi
dx rm -r "dxCompiler_playground:/builds/$username"
dx rm -r "dxCompiler_playground:/builds/user-$username"
dx rm -r "dxCompiler_playground:/unit_tests/$username"
dx rm -r "dxCompiler_playground:/unit_tests/user-$username"

# Run 1 integration test to re-build, upload
if [[ -z "$TEST" ]]; then
  ./scripts/run_tests.py --build only
elif [[ "$TEST" == "--failed" ]]; then
    ./scripts/run_tests.py --failed --delay-compile-errors
else
  ./scripts/run_tests.py --test "$TEST" --delay-compile-errors
fi
