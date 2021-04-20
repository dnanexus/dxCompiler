#!/bin/bash

# Clean artifacts from prev. builds, local
sbt clean && sbt cleanFiles
find . -name target | xargs rm -rf
rm -rf applet_resources
rm dx*.jar

# Clean artifacts from prev. builds, platform
dx rm -r dxCompiler_playground:/builds/`dx whoami`
dx rm -r dxCompiler_playground:/unit_tests/`dx whoami`

# Run 1 integration test to re-build, upload
./scripts/run_tests.py --test add3
