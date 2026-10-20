#!/bin/bash
#
# postProcess.sh AssayName [ProjectName]
#
# Collect the per-pair results under out/ into the project deliverables
# directory <ProjectName>/: the PNGs, the merged IGV .seg (X chromosome
# renamed by fixXChrom.R), the gene calls (getGeneCalls: GeneMatrix.csv,
# SegmentMatrix.csv, <ProjectName>___GeneTable.xlsx) and the PDF of all
# plots (<ProjectName>___seqSeg.pdf). One job in the POST stage of
# PIPE.sh; can be run by hand after selectBestMatch.
#
# The deliverables name each sample by its tumor id only: the tables are
# rewritten by redactDeliverables.R from pairing.tsv (selectBestMatch),
# and the PDF is built from the <pair>_seqSeg.redacted.png plots
# (seqSegment TUMORID=). The PNGs copied to <ProjectName>/ keep the pair
# id; they are for checking, not for delivery.
#
# Every R script runs under the one R 4.5.x on PATH; bin/loadTools.sh
# provides samtools and ImageMagick (magick).
#

set -eu

export SDIR="$( cd "$( dirname "$0" )" && pwd )"

source $SDIR/bin/loadTools.sh

if [ "$#" -lt "1" ]; then
    echo "usage: postProcess.sh AssayName [ProjectName]"
    echo
    echo "Assays"
    ls -1 $SDIR/resources/geneAnnotations | fgrep -v src | awk '{print "   "$1}'
    echo
    exit 1
fi


if [ "$#" -ge "2" ]; then
    projectName=$2
else
    projectName=$($SDIR/extractProjectIDFromPath.py $(pwd -P))
fi
echo
echo projectName=$projectName
echo

assay=$1
if [ ! -e "$SDIR/resources/geneAnnotations/$assay" ]; then
    echo "usage: postProcess.sh AssayName [ProjectName]"
    echo ""
    echo "[$assay] is not a valid assay name"
    echo
    ls -1 $SDIR/resources/geneAnnotations | fgrep -v src | awk '{print "   "$1}'
    echo
    exit 1
fi

PAIRING=pairing.tsv
if [ ! -s "$PAIRING" ]; then
    echo "postProcess.sh: no $PAIRING in $(pwd); run selectBestMatch first"
    exit 1
fi

mkdir -p $projectName
find out -name "*_seqSeg.png" -exec cp {} $projectName \;
find out | fgrep .seg | head -1  | xargs head -1 | cut -f-6 >$projectName/${projectName}___IGV.seg
find out | fgrep .seg | xargs cut -f-6 | fgrep -v "loc.start" >>$projectName/${projectName}___IGV.seg

# Try to infer genome and fix X chromosome
echo
echo "fixXChrom"
echo

Rscript --no-init-file $SDIR/fixXChrom.R $projectName/${projectName}___IGV.seg

echo
echo "getGeneCalls"
echo

find out -type f  | fgrep .seg | perl -pe 's|/[^/]+$|\n|' >lodir

$SDIR/getGeneCalls ASSAY=$assay INPUTS=lodir ODIR=$projectName

echo
echo "redactDeliverables"
echo

Rscript --no-init-file $SDIR/redactDeliverables.R $PAIRING $projectName

echo
echo "seqSeg.pdf"
echo

REDACTED_PNGS=$(find out -name "*_seqSeg.redacted.png" | sort)
if [ -z "$REDACTED_PNGS" ]; then
    echo "postProcess.sh: no *_seqSeg.redacted.png under out/ (seqSegment TUMORID=)"
    exit 1
fi
magick $REDACTED_PNGS $projectName/${projectName}___seqSeg.pdf
