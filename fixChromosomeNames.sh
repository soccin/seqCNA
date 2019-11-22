#!/bin/bash

SDIR="$( cd "$( dirname "$0" )" && pwd )"

IBAM=$1
BASE=$(basename $IBAM | sed 's/.bam//')
OBAM=${BASE}__relabel.bam
ODIR=bamRelabel/${BASE}

mkdir -p $ODIR

if [ -e $ODIR/$OBAM ]; then
    echo "$ODIR/$OBAM already exists, reusing"
    exit
fi

samtools view -H $IBAM | egrep "^@(SQ|HD|RG|PG)" | sed 's/SN:chr/SN:/' >$ODIR/${BASE}.relabel.dict
samtools reheader $ODIR/${BASE}.relabel.dict $IBAM >$ODIR/$OBAM
samtools index $ODIR/$OBAM
