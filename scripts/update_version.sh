#!/usr/bin/env bash

VERSION=$1

# Update the config files with the newest version
CONF_FILES=(
    ./executorWdl/src/main/resources/application.conf
    ./core/src/main/resources/application.conf
    ./executorCwl/src/main/resources/application.conf
    ./compiler/src/main/resources/application.conf
    ./executorCommon/src/main/resources/application.conf
)

# Adapted to work also with macOS version of sed, where -i expects extension arg
# https://stackoverflow.com/questions/19456518/invalid-command-code-despite-escaping-periods-using-sed/19457213#19457213
for i in ${CONF_FILES[@]}; do
    sed -i'' -e "s/version.*$/version = \"${VERSION}\"/" $i
done
