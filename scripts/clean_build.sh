#!/bin/bash

# Clean artifacts from prev. builds, local
sbt clean && sbt cleanFiles
find . -name target | xargs rm -rf
rm -rf applet_resources
rm dx*.jar

# Clean artifacts from prev. builds, platform
username=`dx whoami --id`
dx rm -r dxCompiler_playground:/builds/$username
dx rm -r dxCompiler_playground:/unit_tests/$username

# Run 1 integration test to re-build, upload
./scripts/run_tests.py --test upload_wait
