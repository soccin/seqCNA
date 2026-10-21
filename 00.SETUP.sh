#!/bin/bash
#
# 00.SETUP.sh -- build the repo-local R library once per checkout
#
# Usage, from anywhere:
#    seqCNA/00.SETUP.sh
#
# Builds seqDNAcopy, the fork in the submodule, into Rlib/ (gitignored),
# replacing whatever Rlib/ held. Every other package comes from the
# shared or personal library. R and the packages in Rlib.lock.tsv must
# be at the locked versions; bin/checkRlib.R checks them at the end and
# every R script checks them again (bin/rscript.sh). R itself comes from
# PATH (on IRIS: module load R/4.5.1_nds). Safe to rerun; it rebuilds.
#

set -eu

SDIR=$( cd "$( dirname "$0" )" && pwd )
cd $SDIR

fatal() {
    echo
    echo "    FATAL ERROR: $*"
    echo
    exit 1
}

command -v Rscript >/dev/null 2>&1 || fatal "Rscript not on PATH (module load R/4.5.1_nds)"

LOCK=$SDIR/Rlib.lock.tsv
[ -s "$LOCK" ] || fatal "no $LOCK"

RWANT=$(awk -F'\t' '$1=="R"{print $2}' $LOCK)
RVERSION=$(Rscript --vanilla -e 'cat(as.character(getRversion()))')
if [ "$RVERSION" != "$RWANT" ]; then
    fatal "R $RVERSION on PATH; Rlib.lock.tsv needs R $RWANT (module load R/4.5.1_nds)"
fi
echo "R $RVERSION: $(command -v Rscript)"

if [ ! -e seqDNAcopy/DESCRIPTION ]; then
    fatal "seqDNAcopy/ submodule is empty; run: git submodule update --init"
fi

#
# Compiling goes to node-local scratch, not /tmp (see bin/loadTools.sh)
#
export TMPDIR=/localscratch/$USER/seqCNA_setup.$$
mkdir -p $TMPDIR
trap 'rm -rf $TMPDIR' EXIT

rm -rf Rlib
mkdir Rlib
R CMD INSTALL --clean --preclean -l Rlib seqDNAcopy

#
# Check, and record what was built against
#

export R_LIBS=$SDIR/Rlib${R_LIBS:+:$R_LIBS}

Rscript --vanilla bin/checkRlib.R $LOCK $SDIR/Rlib

{
    Rscript --vanilla bin/checkRlib.R --list $LOCK $SDIR/Rlib
    echo "seqDNAcopy git: $(git -C seqDNAcopy describe --always --dirty)"
    echo "lock md5: $(md5sum <$LOCK | cut -d' ' -f1)"
    echo "built: $(date '+%Y-%m-%d %H:%M:%S %Z') on $(hostname)"
} >$SDIR/Rlib/VERSIONS.txt

cat $SDIR/Rlib/VERSIONS.txt
echo
echo "Rlib/ ready"
