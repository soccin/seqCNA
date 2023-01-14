#!/bin/bash

export SDIR="$( cd "$( dirname "$0" )" && pwd )"

if [ "$#" -lt "3" ]; then
    echo "usage: postProcess.sh AssayName ProjectName WDIR"
    echo
    echo "Assays"
    ls -1 $SDIR/resources/geneAnnotations | fgrep -v src | awk '{print "   "$1}'
    echo
    exit
fi


projectName=$2
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
    exit
fi

WDIR=$3

mkdir -p $projectName

PTAG=$(echo $projectName | tr '/' '_')

find $WDIR/out | fgrep .png | xargs -I % cp % $projectName
find $WDIR/out | fgrep .seg | head -1  | xargs head -1 | cut -f-6 >$projectName/${PTAG}___IGV.seg
find $WDIR/out | fgrep .seg | xargs cut -f-6 | fgrep -v "loc.start" >>$projectName/${PTAG}___IGV.seg

# Try to infer genome and fix X chromosome
echo
echo "fixXChrom"
echo

Rscript --no-save $SDIR/fixXChrom.R $projectName/${PTAG}___IGV.seg

echo
echo "getGeneCalls"
echo

find $WDIR/out -type f  | fgrep .seg | perl -pe 's|/[^/]+$|\n|' >$WDIR/lodir

$SDIR/getGeneCalls ASSAY=$assay INPUTS=$WDIR/lodir ODIR=$projectName

