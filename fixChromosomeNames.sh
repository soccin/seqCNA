#!/bin/bash

set -eu

SDIR="$( cd "$( dirname "$0" )" && pwd )"

module load samtools/1.19.2

IBAM=$1
BASE=$(basename $IBAM | sed 's/.bam//')
OBAM=${BASE}__relabel.bam
ODIR=bamRelabel/${BASE}

mkdir -p $ODIR

if [ -e $ODIR/$OBAM ]; then
    echo "$ODIR/$OBAM already exists, reusing"
    exit
fi

GENOME=$(getGenomeBuildBAM.sh $IBAM)

if [ "$GENOME" != "b37" ]; then

    samtools view -H $IBAM | egrep "^@(SQ|HD|RG|PG)" | sed 's/SN:chr/SN:/' >$ODIR/${BASE}.relabel.dict
    samtools reheader $ODIR/${BASE}.relabel.dict $IBAM >$ODIR/$OBAM
    samtools index $ODIR/$OBAM

else

    #
    # B37 already in the right format just link up
    #

    ln -s $(realpath $IBAM) $ODIR/$OBAM
    if [ -e $(realpath ${IBAM/.bam/.bai}) ]; then
        ln -s $(realpath ${IBAM/.bam/.bai}) $ODIR/${OBAM/.bam/.bai}
    elif [ -e $(realpath ${IBAM}.bai) ]; then
        ln -s $(realpath ${IBAM}.bai) $ODIR/${OBAM/.bam/.bai}
    else
        echo -e "\n\tFATAL ERROR: No BAM INDEX File Found\n"
        echo -e "\t\t$IBAM\n\n"
        exit 1
    fi

fi
