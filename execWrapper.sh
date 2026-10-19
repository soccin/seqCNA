#!/bin/bash

#
# Need to go back to R version 3
#
# First remove current R and R env vars
#
# This is done in PIPE.sh
#
# export PATH=$(echo $PATH | tr ':' '\n' | fgrep -v /R/R | tr '\n' ':')
# unset R_LIBS
# unset R_VERSION

# #
# # Now load R/3.x
# #

# module load R/R-3.6.1
# MAJOROS=7
# export R_VERSION=$(R --version | head -1 | awk '{print $3}')
# export R_LIBS=/home/socci/lib/R/CentOS${MAJOROS}/$R_VERSION

# R --version



export SDIR="$( cd "$( dirname "$0" )" && pwd )"
export SNAME=$(basename $0)

if [ -e /etc/system-release ]; then
    OSSTR=$(cat /etc/system-release)
    OSVER=$(echo $OSSTR | perl -ne 'm/ ((\d+)\.(\d+)(|.\d+)) /; print $1')
    OSMAJOR=$(echo $OSVER | perl -ne 'm/(\d+)\./; print $1')
else
    echo "MISSING /etc/system-release"
    exit 1
fi

case $OSMAJOR in
    6)
    export R_LIBS=$SDIR/Rlib:/home/socci/lib/R/CentOS6/3.4.3
    RSCRIPT=/ifs/work/socci/opt/R/3.4.3/bin/Rscript
    ;;

    7)
    export R_LIBS=$SDIR/Rlib:$R_LIBS
    RSCRIPT=$(which Rscript)
    ;;

    *)
    echo "Unsupportted OS Version" $OSVER, $OSMAJOR
    exit 1
esac

# echo $R_LIBS
# echo $RSCRIPT
# exit

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
