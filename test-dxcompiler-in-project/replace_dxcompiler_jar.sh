#/bin/bash

set -exo pipefail

if [[ `find .. -maxdepth 1 -name "dxCompiler-*.jar" | wc -l` -eq 1 ]]; then
    echo "Copying dxCompiler jar from root directory to applet resources."
    cp ../dxCompiler-*.jar resources/home/dnanexus/dxCompiler.jar
else
    echo "Ensure that 1 copy of the dxCompiler jar is staged in the root directory."
    exit 1
fi
