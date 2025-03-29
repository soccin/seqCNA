#!/bin/bash

if [ "$#" != "1" ]; then
    echo
    echo "    usage: delivery.sh /ifs/res/seq/pi/invest/r_###"
    echo
    exit
fi

PDIR=$(realpath $1)

#echo $PDIR

BASE=$(basename $PDIR)
PNAME=$(basename $(dirname $PDIR))
echo $PDIR, $PNAME

if [[ ! $BASE =~ r_[0-9]+ ]] || [[ ! $PNAME =~ Proj_ ]]; then
    echo "Invalid propect delivery directory"
fi

mkdir $PDIR/seqcna

cp $PNAME/*csv $PDIR/seqcna
cp $PNAME/*___seqSeg.pdf $PDIR/seqcna
cp $PNAME/*___IGV.seg $PDIR/seqcna
cp $PNAME/P*___GeneTable.xlsx $PDIR/seqcna
