#!/bin/bash
#
# execWrapper.sh -- run the R script that shares the invoking name
#
# getPairedCounts, seqSegment, getGeneCalls and selectBestMatch are
# symlinks to this file. It exports SDIR (this checkout) and SNAME (the
# invoked name), puts the repo-local Rlib/ (built by 00.SETUP.sh) first
# on R_LIBS, and runs Rscript --vanilla $0.R. With no arguments it
# prints docs/$SNAME.doc instead.
#
# R 4.5.x is taken from PATH (on IRIS: module load R/4.5.1_nds). The
# common packages come from the shared or personal library; --vanilla
# still honors the R_LIBS and R_LIBS_USER environment variables, only
# ~/.Rprofile and ~/.Renviron are skipped.
#

export SDIR="$( cd "$( dirname "$0" )" && pwd )"
export SNAME=$(basename $0)

fatal() {
    echo
    echo "    FATAL ERROR: $SNAME: $*"
    echo
    exit 1
}

RSCRIPT=$(command -v Rscript) || fatal "Rscript not on PATH (module load R/4.5.1_nds)"

case $($RSCRIPT --vanilla -e 'cat(as.character(getRversion()))') in
    4.5.*) ;;
    *) fatal "need R 4.5.x on PATH, found $($RSCRIPT --version 2>&1 | head -1)" ;;
esac

if [ ! -e $SDIR/Rlib/seqDNAcopy/DESCRIPTION ]; then
    fatal "$SDIR/Rlib/seqDNAcopy missing; run $SDIR/00.SETUP.sh"
fi

export R_LIBS=$SDIR/Rlib${R_LIBS:+:$R_LIBS}

DOCFILE=$SDIR/docs/${SNAME}.doc
if [ "$#" == "0" ] && [ -e $DOCFILE ]; then
    # If no args then just print usage
    # do not bother running R script

    cat $DOCFILE
    echo
    echo
    exit

fi

exec $RSCRIPT --vanilla --no-save "$0.R" "$@"
