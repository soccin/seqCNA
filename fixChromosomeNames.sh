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

ln -s $(realpath $IBAM) $ODIR
ln -s $(realpath ${IBAM/.bam/.bai}) $ODIR
