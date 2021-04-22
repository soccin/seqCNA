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
    echo "Invalid project delivery directory"
fi

ODIR=$(basename $(dirname $(find . | fgrep ___GeneTable.xlsx)))

if [ "$PNAME" == "$ODIR" ]; then

    mkdir -p $PDIR/seqcna

    cp $ODIR/*csv $PDIR/seqcna
    cp $ODIR/*___seqSeg.pdf $PDIR/seqcna
    cp $ODIR/*___IGV.seg $PDIR/seqcna
    cp $ODIR/P*___GeneTable.xlsx $PDIR/seqcna

else

    mkdir -p $PDIR/seqcna/$ODIR

    cp $ODIR/*csv $PDIR/seqcna/$ODIR
    cp $ODIR/*___seqSeg.pdf $PDIR/seqcna/$ODIR
    cp $ODIR/*___IGV.seg $PDIR/seqcna/$ODIR
    cp $ODIR/P*___GeneTable.xlsx $PDIR/seqcna/$ODIR

fi

