#!/bin/bash

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

    8)
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
