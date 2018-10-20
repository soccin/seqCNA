#!/bin/bash

samtools view -H $1 \
    | egrep "^@RG" \
    | tr '\t' '\n' \
    | fgrep SM: \
    | sed 's/SM://' 

