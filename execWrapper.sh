#!/bin/bash
#
# execWrapper.sh -- run the R script that shares the invoking name
#
# getPairedCounts, seqSegment, getGeneCalls and selectBestMatch are
# symlinks to this file. It exports SDIR (this checkout) and SNAME (the
# invoked name) and runs $0.R through bin/rscript.sh, which puts Rlib/
# (seqDNAcopy, built by 00.SETUP.sh) first on R_LIBS and checks R and
# the packages in Rlib.lock.tsv first. With no arguments it prints
# docs/usage/$SNAME.doc instead.
#

export SDIR="$( cd "$( dirname "$0" )" && pwd )"
export SNAME=$(basename $0)

DOCFILE=$SDIR/docs/usage/${SNAME}.doc
if [ "$#" == "0" ] && [ -e $DOCFILE ]; then
    # If no args then just print usage
    # do not bother running R script

    cat $DOCFILE
    echo
    echo
    exit

fi

exec $SDIR/bin/rscript.sh "$0.R" "$@"
