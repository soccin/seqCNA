#!/bin/bash

SID=$(samtools view -H $1 \
    | egrep "^@RG" \
    | tr '\t' '\n' \
    | fgrep SM: \
    | sed 's/SM://' \
    | uniq \
    | xargs \
    | tr ' ' ',')

if [[ ! $SID =~ ^s_ ]]; then
    echo "s_"$SID
else
    echo $SID
fi

