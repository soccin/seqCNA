#!/bin/bash

normal=$1
tumor=$2

tumorId=$(basename $tumor | sed 's/.bam//')

if [ "$#" == "3" ]; then
    sID=$3
else
    sID=${tumorId}__$(basename $normal | sed 's/.bam//')
fi

scatter=$(echo $tumorId | perl -ne 'print substr($_,0,2)')
mkdir -p LSF/$scatter

oDir=out/$scatter/$tumorId/$sID

bsub -m commonHG -o LSF/$scatter -J WGSCNA_$sID -We 59 -R "rusage[iounits=.1]" \
    ./seqCNA/getPairedCounts NORMAL=$normal TUMOR=$tumor \
        ODIR=$oDir \
        SAMPLEID=$sID

bsub -m commonHG -o LSF/$scatter -J SEQSEG_$sID -We 59 -R "rusage[iounits=.1]" -w "post_done(WGSCNA_$sID)" \
    ./seqCNA/seqSegment BINSIZE=100 \
        ODIR=$oDir \
        COUNTS=$oDir/${sID}_Counts.rda

#bsub -m commonHG -o LSF/$scatter -J CLEANUP_$sID -We 59 -w "post_done(SEQSEG_$sID)" \
#    rm $oDir/${sID}_Counts.rda $oDir/${sID}_seqSeg.rda

