#!/bin/bash

if [ "$#" -lt "2" ]; then
    echo "usage: seqCNA.sh BINSIZE normal.bam tumor.bam [SampleID]"
    exit
fi

SDIR="$( cd "$( dirname "$0" )" && pwd )"

BINSIZE=$1
normal=$2
tumor=$3

tumorId=$(basename $tumor | sed 's/.bam//')

GENOME_TAG=$($SDIR/GenomeData/getGenomeBuildBAM.sh $tumor)

case $GENOME_TAG in
    mm10-relabel)
    GENOME="mm10"
    ;;

    GRC_m38)
    GENOME="mm10"
    ;;

    b37)
    GENOME="hg19"
    ;;

    *)
    echo "Unknown or un-usable genome"
    echo $GENOME_TAG
    exit
    ;;
esac

if [ "$#" == "4" ]; then
    sID=$4
else
    sID=${tumorId}__$(basename $normal | sed 's/.bam//')
fi

scatter=$(echo $tumorId | perl -ne 'print substr($_,0,2)')
mkdir -p LSF/$scatter

oDir=out/$scatter/$tumorId/$sID

bsub -m commonHG -o LSF/$scatter -J WGSCNA_$sID -We 119 -R "rusage[mem=32]" \
    $SDIR/getPairedCounts GENOME=$GENOME NORMAL=$normal TUMOR=$tumor \
        ODIR=$oDir \
        SAMPLEID=$sID

bsub -m commonHG -o LSF/$scatter -J SEQSEG_$sID -We 119 -R "rusage[mem=32]" -w "post_done(WGSCNA_$sID)" \
    $SDIR/seqSegment \
        BINSIZE=$BINSIZE \
        ODIR=$oDir \
        COUNTS=$oDir/${sID}_Counts.rda

#bsub -m commonHG -o LSF/$scatter -J CLEANUP_$sID -We 59 -w "post_done(SEQSEG_$sID)" \
#    rm $oDir/${sID}_Counts.rda $oDir/${sID}_seqSeg.rda

