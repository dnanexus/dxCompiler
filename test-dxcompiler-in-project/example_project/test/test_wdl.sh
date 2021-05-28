#/bin/bash

set -exo pipefail

DXCOMPILER_JAR=/home/dnanexus/dxCompiler.jar

# Download files needed for compiling workflows
dx download -r workflow_1

# Set up test output folder for each workflow
dx mkdir -p test_out/workflow_1

# Compile workflow 1
exec_1=`java -jar $DXCOMPILER_JAR compile workflow_1/workflow_1.wdl \
	-folder /test_out/workflow_1 \
	-force`

# Run workflow 1 if compiled
# Use detached jobs so workflows fail individually
if [[ "$exec_1" == "workflow-"* ]]; then
    dx run $exec_1 \
        --destination /test_out/workflow_1 \
        --detach \
        -y
fi
