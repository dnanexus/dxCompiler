#!/bin/bash

PROJECT_ID_FILE=$1
INVITEE=$2
LEVEL=$3

for i in `cat $PROJECT_ID_FILE | awk '{ print $2}'`; do
    echo $i
    dx api $i invite '{"invitee": "'$INVITEE'", "level": "'$LEVEL'"}'
done
