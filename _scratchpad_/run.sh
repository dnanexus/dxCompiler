#!/bin/bash

PROJECT=project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq
FOLDER=/rkepych
DXCOMPILER_JAR=../dxCompiler-2.11.6-SNAPSHOT.jar

RUN_OPTS='--ignore-reuse --destination ${PROJECT}:${FOLDER} --preserve-job-'

workflow=$(java -jar $DXCOMPILER_JAR compile prepare_scatter_gather.wdl -project $PROJECT -folder $FOLDER -force)

dx run $workflow $RUN_OPTS
