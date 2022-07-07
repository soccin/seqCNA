QSYNC=#


##
# QRUN NCORES QTAG <HOLD hold_id> <VMEM size> LONG
#
# HOLD, VMEM and LONG are optional but if you more than one
# they must be in this order
#   HOLD
#   VMEM
#   LONG
#

QRUN () {

    NCORES=$1
    QTAG=$2
    echo QTAG=$QTAG
    shift 2
    QHOLD=""
    if [ "$1" == "HOLD" ]; then
        QHOLD="-w post_done($2)"
        shift 2
        echo QHOLD=$QHOLD
    fi

    VMEM=""
    if [ "$1" == "VMEM" ]; then

        TOTALMEM=$2
        MEMPERSLOT=$((TOTALMEM / NCORES))
        VMEM='-R "rusage[mem='$MEMPERSLOT']"'

        shift 2
        echo VMEM=$VMEM

    fi

    TIME=$TIME_SHORT
    if [ "$1" == "LONG" ]; then
        TIME=$TIME_LONG
        shift 1
        echo LONG Job
    fi

    #RET=$(bsub $TIME $QHOLD $VMEM -n $NCORES -J $QTAG -o LSF.PEMAP/ $*)
    echo RET=bsub $QHOLD $VMEM -n $NCORES -J $QTAG -o LSF.PEMAP/ $*
    echo "#QRUN RET=" $RET
    echo

}
