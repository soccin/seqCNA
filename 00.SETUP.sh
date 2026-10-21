#!/bin/bash
#
# 00.SETUP.sh -- build the repo-local R library once per checkout
#
# Usage, from anywhere:
#    seqCNA/00.SETUP.sh
#
# Installs every R package the pipeline uses, at the versions in
# Rlib.lock.tsv, into Rlib/ (gitignored): each from its source tarball
# (CRAN, Bioconductor or GitHub, md5-checked), then seqDNAcopy from the
# submodule. Nothing is taken from the user or site library. R itself
# comes from PATH and must be the locked version (on IRIS: module load
# R/4.5.1_nds).
#
# The library is built in Rlib.new/ and replaces Rlib/ only when it is
# complete and passes bin/checkRlib.R, so a failed build leaves the old
# one in place. The tarballs are kept in Rlib/src/ and reused by the
# next build. Safe to rerun; it rebuilds. Takes about half an hour,
# mostly compiling. SETUP_NCPU (default 4) sets make -j.
#
# To change a version, write a new lock with bin/mkRlibLock.R (see
# there); never install into Rlib/ by hand.
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

RWANT=$(awk -F'\t' '$3=="R"{print $2}' $LOCK)
RVERSION=$(Rscript --vanilla -e 'cat(as.character(getRversion()))')
if [ "$RVERSION" != "$RWANT" ]; then
    fatal "R $RVERSION on PATH; Rlib.lock.tsv needs R $RWANT (module load R/4.5.1_nds)"
fi
echo "R $RVERSION: $(command -v Rscript)"

if [ ! -e seqDNAcopy/DESCRIPTION ]; then
    fatal "seqDNAcopy/ submodule is empty; run: git submodule update --init"
fi

STAGE=$SDIR/Rlib.new
rm -rf $STAGE
mkdir -p $STAGE/src

#
# Compiling goes to node-local scratch, not /tmp (see bin/loadTools.sh)
#
export TMPDIR=/localscratch/$USER/seqCNA_setup.$$
mkdir -p $TMPDIR
trap 'rm -rf $TMPDIR' EXIT

export R_LIBS=$STAGE
export R_LIBS_USER=$STAGE
export R_LIBS_SITE=$STAGE
export MAKEFLAGS=-j${SETUP_NCPU:-4}

#
# Locked packages, then seqDNAcopy from the submodule
#

Rscript --vanilla bin/installRlib.R $LOCK $STAGE $STAGE/src $SDIR/Rlib/src

R CMD INSTALL --clean --preclean -l $STAGE seqDNAcopy

Rscript --vanilla bin/checkRlib.R $LOCK $STAGE

rm -rf $SDIR/Rlib.old
if [ -d $SDIR/Rlib ]; then
    mv $SDIR/Rlib $SDIR/Rlib.old
fi
mv $STAGE $SDIR/Rlib
rm -rf $SDIR/Rlib.old

#
# Record what was built
#

export R_LIBS=$SDIR/Rlib
export R_LIBS_USER=$SDIR/Rlib
export R_LIBS_SITE=$SDIR/Rlib

{
    Rscript --vanilla bin/checkRlib.R --list $LOCK $SDIR/Rlib
    echo "seqDNAcopy git: $(git -C seqDNAcopy describe --always --dirty)"
    echo "lock md5: $(md5sum <$LOCK | cut -d' ' -f1)"
    echo "built: $(date '+%Y-%m-%d %H:%M:%S %Z') on $(hostname)"
} >$SDIR/Rlib/VERSIONS.txt

cat $SDIR/Rlib/VERSIONS.txt
echo
echo "Rlib/ ready"
