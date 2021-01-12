#!/bin/bash

PROJECT_ID_FILE=$1

for i in `cat $PROJECT_ID_FILE | awk '{ print $2}'`; do
    echo $i
    dx api $i invite '{"invitee": "PUBLIC", "level": "VIEW"}'
done
