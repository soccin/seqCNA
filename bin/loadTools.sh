#!/bin/bash
#
# loadTools.sh -- put the external tools the pipeline needs on PATH
#
# Usage: . $SDIR/bin/loadTools.sh
#
# samtools comes from Lmod. Batch jobs inherit the exported module
# function through --export=ALL, but a bare shell may not have it, so
# source the init script when it is missing.
#
# R 4.5.x is expected on PATH already (on IRIS: module load
# R/4.5.1_nds, a personal module in ~/.modules); the wrapped R tools
# check it themselves, this only reports it early.
#
# ImageMagick (magick, for the PNG -> PDF merge in postProcess.sh) has
# no system module; it is prepended from the core001 shared tree.
#
# Safe to source under set -u.
#

SEQCNA_SAMTOOLS_MODULE=${SEQCNA_SAMTOOLS_MODULE:-samtools/1.20}
SEQCNA_IMAGEMAGICK_BIN=${SEQCNA_IMAGEMAGICK_BIN:-/usersoftware/core001/common/RHEL_8/imagemagick/7.1.2-33/bin}

if [ "$(type -t module)" != "function" ]; then
    . /etc/profile.d/modules.sh
fi

module load $SEQCNA_SAMTOOLS_MODULE

if [[ ":$PATH:" != *":$SEQCNA_IMAGEMAGICK_BIN:"* ]]; then
    export PATH=$SEQCNA_IMAGEMAGICK_BIN:$PATH
fi

for SEQCNA_TOOL in samtools magick Rscript; do
    if ! command -v $SEQCNA_TOOL >/dev/null 2>&1; then
        echo
        echo "    FATAL ERROR: $SEQCNA_TOOL not on PATH"
        echo
        exit 1
    fi
done
unset SEQCNA_TOOL

case $(Rscript --vanilla -e 'cat(as.character(getRversion()))') in
    4.5.*) ;;
    *)
        echo
        echo "    FATAL ERROR: need R 4.5.x on PATH (module load R/4.5.1_nds)"
        echo "    found: $(command -v Rscript)"
        echo
        exit 1
        ;;
esac

#
# Temp files never go to /tmp. IRIS compute nodes set TMPDIR=/tmp, a
# 137G volume shared by every job on the node; the R session tempdir
# and samtools sort spill land there by default. Point TMPDIR at the
# node-local /localscratch (2.8T) instead. bin/runStage.sh has already
# set a per-job directory under it for stage jobs; keep that. This also
# covers tools run by hand.
#

SEQCNA_LOCAL_TMP=${SEQCNA_LOCAL_TMP:-/localscratch/$USER}

if [[ "${TMPDIR:-}" != "$SEQCNA_LOCAL_TMP"/* ]]; then
    export TMPDIR=$SEQCNA_LOCAL_TMP
fi

if ! mkdir -p "$TMPDIR"; then
    echo
    echo "    FATAL ERROR: cannot create TMPDIR=[$TMPDIR]"
    echo
    exit 1
fi
