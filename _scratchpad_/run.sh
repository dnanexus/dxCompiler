#!/bin/bash

PROJECT=project-Fy9QqgQ0yzZbg9KXKP4Jz6Yq
FOLDER=/rkepych

java -jar dxCompiler.jar compile prepare_scatter_gather.wdl -project $PROJECT -folder $FOLDER
