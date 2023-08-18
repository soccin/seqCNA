#!/bin/bash

export SDIR="$( cd "$( dirname "$0" )" && pwd )"

if [ "$#" -lt "1" ]; then
    echo "usage: postProcess.sh AssayName [ProjectName]"
    echo
    echo "Assays"
    ls -1 $SDIR/resources/geneAnnotations | fgrep -v src | awk '{print "   "$1}'
    echo
    exit
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
    exit
fi

mkdir $projectName
find out | fgrep .png | xargs -I % cp % $projectName
find out | fgrep .seg | head -1  | xargs head -1 | cut -f-6 >$projectName/${projectName}___IGV.seg
find out | fgrep .seg | xargs cut -f-6 | fgrep -v "loc.start" >>$projectName/${projectName}___IGV.seg

# Try to infer genome and fix X chromosome
echo
echo "fixXChrom"
echo

Rscript_4x=/juno/res/bic/shared/Multiomyx/Projects/GBM/opt/bin/Rscript
$Rscript_4x $SDIR/fixXChrom.R $projectName/${projectName}___IGV.seg

echo
echo "getGeneCalls"
echo

find out -type f  | fgrep .seg | perl -pe 's|/[^/]+$|\n|' >lodir

$SDIR/getGeneCalls ASSAY=$assay INPUTS=lodir ODIR=$projectName

