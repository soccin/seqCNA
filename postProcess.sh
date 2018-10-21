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
    projectName=$(basename $PWD)
fi

assay=$1

mkdir $projectName
find out | fgrep .png | xargs -I % cp % $projectName
find out | fgrep .seg | head -1  | xargs head -1 | cut -f-6 >$projectName/${projectName}___IGV.seg
find out | fgrep .seg | xargs cut -f-6 | fgrep -v "loc.start" >>$projectName/${projectName}___IGV.seg

find out -type f  | fgrep .seg | perl -pe 's|/[^/]+$|\n|' >lodir
$SDIR/getGeneCalls ASSAY=$assay INPUTS=lodir ODIR=$projectName

