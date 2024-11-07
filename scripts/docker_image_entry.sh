#!/bin/sh
echo "WARNING: This Docker image will be removed permanently on 2024/12/1. Please STOP using this image ASAP"
exec java -jar /dxCompiler.jar "$@"