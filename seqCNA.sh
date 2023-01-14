#!/bin/bash

if [ "$#" -lt "4" ]; then
    echo "usage: seqCNA.sh BINSIZE normal.bam tumor.bam WDIR"
    echo "NUM=" $#
    exit
fi

SDIR="$( cd "$( dirname "$0" )" && pwd )"

BINSIZE=$1
normal=$2
tumor=$3
WDIR=$4

tumorId=$($SDIR/getSampleTag.sh $tumor)

GENOME_TAG=$($SDIR/GenomeData/getGenomeBuildBAM.sh $tumor)

case $GENOME_TAG in
    mm10-relabel)
    GENOME="mm10"
    ;;

    GRC_m38)
    GENOME="mm10"
    ;;

    b37+mm10)
    GENOME="hg19"
    ;;

    b37)
    GENOME="hg19"
    ;;

    b37_dmp)
    GENOME="hg19"
    ;;

    *)
    echo "Unknown or un-usable genome"
    echo $GENOME_TAG
    exit 13
    ;;
esac

scatter=$(echo $tumorId | md5sum - | perl -ne 'print substr($_,0,2)')
LSFDIR=LSF/$scatter/$tumorId
mkdir -p $LSFDIR

sID=${tumorId}__$($SDIR/getSampleTag.sh $normal)

oDir=$WDIR/out/$scatter/$tumorId/$sID

bsub -o $LSFDIR -J ${QTAG}_WGSCNA_$sID -W 59 -n 2 -R "rusage[mem=16]" \
    $SDIR/getPairedCounts GENOME=$GENOME NORMAL=$normal TUMOR=$tumor \
        ODIR=$oDir \
        SAMPLEID=$sID

bsub -o $LSFDIR -J ${QTAG}_SEQSEG_$sID -W 59 -n 2 -R "rusage[mem=16]" -w "post_done(${QTAG}_WGSCNA_$sID)" \
    $SDIR/seqSegment \
        BINSIZE=$BINSIZE \
        ODIR=$oDir \
        COUNTS=$oDir/${sID}_Counts.rda

#bsub -o LSF/$scatter -J CLEANUP_$sID -W 59 -w "post_done(SEQSEG_$sID)" \
#    rm $oDir/${sID}_Counts.rda $oDir/${sID}_seqSeg.rda

