#!/bin/bash

OROOT="out"

getInsStats () {

    bi=$1

    md5=$(echo $bi | md5sum | awk '{print $1}')

    D1=$(echo $md5 | perl -ne 'm/(..)/; print $1')
    D2=$(echo $md5 | perl -ne 'm/..(..)/; print $1')

    BASE=$(basename ${bi/.bam/})

    ODIR=$OROOT/$D1/$D2/$BASE
    mkdir -p $ODIR

    picardV2 CollectInsertSizeMetrics \
        I=$bi H=$ODIR/${BASE}___INS.pdf O=$ODIR/${BASE}___INS.txt

}

for bam in $*; do

    getInsStats $bam

done
