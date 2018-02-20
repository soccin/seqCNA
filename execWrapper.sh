#!/bin/bash

export SDIR="$( cd "$( dirname "$0" )" && pwd )"
export SNAME=$(basename $0)
export R_LIBS=$SDIR/Rlib:/home/socci/lib/R/3.3.1
RSCRIPT=/opt/common/CentOS_6-dev/R/R-3.3.1/bin/Rscript

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
