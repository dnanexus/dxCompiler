#/bin/bash

set -exo pipefail

DXCOMPILER_JAR=/home/dnanexus/dxCompiler.jar

# Timestamp for unique test run folder names
TIMESTAMP=$(date +%Y%m%d%H%M%S)

# Download files needed for compiling workflows
dx download -r workflow_1

# Set up test output folder for each workflow
WF1_DIR=/test_out/workflow_1_$TIMESTAMP
dx mkdir -p $WF1_DIR

# Compile workflow 1
exec_1=`java -jar $DXCOMPILER_JAR compile workflow_1/workflow_1.wdl \
    -extras workflow_1/extras.json \
    -folder $WF1_DIR \
    -force`

# Run workflow 1 if compiled
# Use detached jobs so workflows fail individually
# Ignore reuse so full workflow runs each time, if benchmarking
if [[ "$exec_1" == "workflow-"* ]]; then
    dx run $exec_1 \
        --destination $WF1_DIR \
        --detach \
        --ignore-reuse \
        --brief
fi
