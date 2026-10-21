#!/bin/bash
#
# rscript.sh -- run an R script on the pinned library
#
# Usage:
#    seqCNA/bin/rscript.sh SCRIPT.R [ARGS...]
#    seqCNA/bin/rscript.sh --check
#
# Every R script of the pipeline runs through this: execWrapper.sh for
# the wrapped tools, postProcess.sh for fixXChrom.R and
# redactDeliverables.R, and by hand for plotSeg.R. R's library path is
# set to Rlib/ (built by 00.SETUP.sh) and R's own library only, so
# nothing is taken from a user or site library, and bin/checkRlib.R
# checks R and every package against Rlib.lock.tsv before the script
# runs; a mismatch stops it with exit 1. The script runs under
# Rscript --vanilla (no ~/.Rprofile, no ~/.Renviron).
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

#
# R_HOME/etc/Renviron sets each of these only when it is unset or
# empty, so all three must be set, to Rlib/.
#
export R_LIBS=$RS_SDIR/Rlib
export R_LIBS_USER=$RS_SDIR/Rlib
export R_LIBS_SITE=$RS_SDIR/Rlib

$RSCRIPT --vanilla $RS_SDIR/bin/checkRlib.R $RS_SDIR/Rlib.lock.tsv $RS_SDIR/Rlib || exit 1

if [ "$1" == "--check" ]; then
    exit 0
fi

exec $RSCRIPT --vanilla "$@"
