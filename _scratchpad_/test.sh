#!/bin/bash

i=0
while [[ $i -le 5 ]]; do
    echo $((2**i))
    sleep $((2**i))
    i=$((i + 1))
done
