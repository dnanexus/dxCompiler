#!/bin/bash
rm -Rf corpora
mkdir corpora
jq '.corpora | .[] | "cd corpora && repo=\(.url | split("/")[4] | split(".")[0]) && if [ ! -d $repo ] ; then git clone --recurse-submodules \(.url) && cd $repo && git checkout \(.tag) ; fi"' corpora_repos.json | xargs -i -t bash -c "{}"