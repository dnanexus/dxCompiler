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
for i in ${CONF_FILES[@]}; do
    sed -i'' -e "s/version.*$/version = \"${VERSION}\"/" $i
done