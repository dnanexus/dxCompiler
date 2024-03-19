#!/bin/bash

PROJECT=project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq
FOLDER=/rkepych
DXCOMPILER_JAR=~/Documents/dxCompiler/applet_resources/dxCompiler.jar

java -jar $DXCOMPILER_JAR compile prepare_scatter_gather.wdl -project $PROJECT -folder $FOLDER
