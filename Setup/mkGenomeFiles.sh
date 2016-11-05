#!/bin/bash

SDIR="$( cd "$( dirname "$0" )" && pwd )"

GENOME=$1
BINSIZE=$2
GENOMEFASTA=$3

GTAG=${GENOME/.genome/}

bedtools makewindows -g $GENOME -w $BINSIZE >${GTAG}.BINS_${BINSIZE}.bed

bedtools makewindows -w 1000 -s 100 -g $GENOME \
    | bedtools nuc -fi $GENOMEFASTA -bed - \
    | gzip -c - > ${GTAG}_1000by100_TileBin.nuc.gz

Rscript --no-save $SDIR/mkRDA.R ${GTAG}_1000by100_TileBin.nuc.gz

