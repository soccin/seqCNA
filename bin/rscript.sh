#!/bin/bash
#
# rscript.sh -- run an R script after checking the locked versions
#
# Usage:
#    seqCNA/bin/rscript.sh SCRIPT.R [ARGS...]
#    seqCNA/bin/rscript.sh --check
#
# Every R script of the pipeline runs through this: execWrapper.sh for
# the wrapped tools, postProcess.sh for fixXChrom.R and
# redactDeliverables.R, and by hand for plotSeg.R. Rlib/ (seqDNAcopy,
# built by 00.SETUP.sh) goes first on R_LIBS; the other packages come
# from the shared or personal library. bin/checkRlib.R checks R and the
# packages in Rlib.lock.tsv before the script runs; a mismatch stops it
# with exit 1. The script runs under Rscript --vanilla (no ~/.Rprofile,
# no ~/.Renviron; R_LIBS and R_LIBS_USER from the environment still
# apply).
#
# --check runs only the check.
#
# R itself is taken from PATH (on IRIS: module load R/4.5.1_nds).
#

RS_SDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"

RSCRIPT=$(command -v Rscript) || {
    echo
    echo "    FATAL ERROR: Rscript not on PATH (module load R/4.5.1_nds)"
    echo
    exit 1
}

export R_LIBS=$RS_SDIR/Rlib${R_LIBS:+:$R_LIBS}

$RSCRIPT --vanilla $RS_SDIR/bin/checkRlib.R $RS_SDIR/Rlib.lock.tsv $RS_SDIR/Rlib || exit 1

if [ "$1" == "--check" ]; then
    exit 0
fi

exec $RSCRIPT --vanilla "$@"
