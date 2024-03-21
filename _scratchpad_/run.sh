#!/bin/bash

PROJECT=project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq
FOLDER=/rkepych
DXCOMPILER_JAR=../dxCompiler-2.11.6-SNAPSHOT.jar
WDL_SRC=scatter_collect_with_glob.wdl

RUN_OPTS="--ignore-reuse --destination ${PROJECT}:${FOLDER}"

workflow=$(java -jar $DXCOMPILER_JAR compile $WDL_SRC -project $PROJECT -folder $FOLDER -force)

dx run $workflow $RUN_OPTS
