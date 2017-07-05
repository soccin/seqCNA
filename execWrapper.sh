#!/bin/bash

export SDIR="$( cd "$( dirname "$0" )" && pwd )"
export R_LIBS=$SDIR/Rlib:$R_LIBS

DOCFILE=$SDIR/docs/$(basename $0).doc
if [ "$#" == "0" ] && [ -e $DOCFILE ]; then
    # If no args then just print usage
    # do not bother running R script

    cat $DOCFILE
    echo
    echo
    exit

fi

exec /opt/common/CentOS_6-dev/R/R-3.2.2/bin/Rscript --vanilla --no-save "$0.R" "$@"
