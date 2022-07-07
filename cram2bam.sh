#!/bin/bash

SDIR="$( cd "$( dirname "$0" )" && pwd )"

ICRAM=$1
BASE=$(basename $ICRAM | sed 's/.cram//')
OBAM=${BASE}
ODIR=cvtBams/${BASE}

mkdir -p $ODIR

if [ -e $ODIR/$OBAM ]; then
    echo "$ODIR/$OBAM already exists, reusing"
    exit
fi

BUILD=$($SDIR/GenomeData/getGenomeBuildBAM.sh $ICRAM)

if [ -e $SDIR/GenomeData/genomeInfo_${BUILD}.sh ]; then
    . $SDIR/GenomeData/genomeInfo_${BUILD}.sh
else
    echo
    echo "Invalid Genome Build" $BUILD
    echo
    exit
fi

samtools view -T $GENOME -o $ODIR/$BAM/$OBAM $ICRAM
#samtools index $ODIR/$OBAM
