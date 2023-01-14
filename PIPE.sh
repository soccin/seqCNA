#!/bin/bash

#git clone git@github.com:soccin/seqCNA.git

#
# Need to go back to R version 3
#
# First remove current R and R env vars

export PATH=$(echo $PATH | tr ':' '\n' | fgrep -v /R/R | tr '\n' ':')
unset R_LIBS
unset R_VERSION

#
# Now load R/3.x
#

module load R/R-3.6.1
MAJOROS=7
export R_VERSION=$(R --version | head -1 | awk '{print $3}')
export R_LIBS=/home/socci/lib/R/CentOS${MAJOROS}/$R_VERSION

R --version

SDIR="$( cd "$( dirname "$0" )" && pwd )"

function usage {
    echo
    echo "   usage: $(dirname $SDIR)/PIPE.sh [sWGS|TARGETTED] DB_OF_NORMALS.csv DB_NORMAL_DECILES.csv TARGET_BAM"
    echo "       sWGS                Shallow Whole Genome Seq"
    echo "       TARGETTED           Targeted assays (IMPACT/EXOME/...)"
    echo "       DB_OF_NORMALS       CSV file or NORMAL BAMS"
    echo "       DB_NORMAL_DECILES   CSV file of NORMAL DECILES"
    echo "       TARGET_BAM          BAM File Path of Target (Tumor)"
    echo
    exit
}

if [ "$#" -ne 4 ]; then
    usage
fi

TAG=qSeqCNA
UUID=$(uuidgen -t | cut -d- -f-4)
export QTAG=${TAG}_${UUID}

ASSAY=$(echo $1 | tr '[a-z]' '[A-Z]')
if [ "$ASSAY" == "SWGS" ]; then
    BINSIZE=auto
    echo
    echo "Setting BINSIZE=auto"
elif [ "$ASSAY" == "TARGETTED" ]; then
    BINSIZE=100
    echo
    echo "Setting BINSIZE=100"
else
    echo
    echo "Unknown assay type ["$ASSAY"]"
    usage
fi

NORMAL_DB=$2
DECILE_DB=$3
TUMOR=$4

SID=$($SDIR/getSampleTag.sh $TUMOR)

DTS=$(date +%Y%m%d%H%M%S)
WDIR=_scratch/$DTS/$SID
mkdir -vp $WDIR

echo "=============================================================================="
echo "Finding Best Decile Normals"
echo
echo

$SDIR/getBestDecileNormals $DECILE_DB $TUMOR $WDIR

fgrep -wf <(cat $WDIR/bestNormals___$SID | awk '{print $1","}') $NORMAL_DB | cut -f2 -d, >$WDIR/normalBams

echo "=============================================================================="
echo "Running $SDIR/seqCNA.sh"
echo
echo

for normal in $(cat $WDIR/normalBams); do
    $SDIR/seqCNA.sh $BINSIZE $normal $TUMOR
    EXITCODE=$?
    if [ "$EXITCODE" != "0" ]; then
        echo
        echo "ERROR in "$SDIR/seqCNA.sh $BINSIZE $normal $tumor
        echo "CODE = "$EXITCODE
        echo
        exit 1
    fi
done

exit

bSync "${QTAG}_.*"

ERR2=$(parseLSF.py LSF/*/* | fgrep -v Succ)
if [ "$ERR2" != "" ]; then
    echo "ERROR @ Stage2"
    parseLSF.py LSF/*/*| fgrep -v Succ
    exit
fi

echo "=============================================================================="
echo "Running $SDIR/selectBestMatch"
echo
echo

$SDIR/selectBestMatch out
# Take best match from T/N Pairs
#cat pipeline/*_sample_pairing.txt | awk '{print $2"__"$1}' | fgrep -v POOL >bestMatches____out

mkdir outAll
rsync -avP --link-dest=../out out/ outAll

ls -d out/s_/*/* | fgrep -vf bestMatches____out | xargs -t rm -rf

GENOME=$($SDIR/GenomeData/getGenomeBuildBAM.sh $(ls $BAMDIR/*.bam | head -1))
find bamRelabel | fgrep .bam | head -1
if [ -e pipeline/*request.txt ]; then
    PROJNO=$(echo $(ls pipeline/*request.txt) | perl -ne 'm|/(Proj_.*)_request|; print $1')
else
    PROJNO=$($SDIR/extractProjectIDFromPath.py $(pwd -P))
fi

if [[ $GENOME =~ b37|hg19 ]]; then
    ASSAY=Exome
elif [[ $GENOME =~ ^mm10 ]]; then
    ASSAY=M-IMPACT_v1
else
    echo
    echo "Unknown geneome"
    echo "   "$GENOME
    echo
    exit -1
fi

echo "=============================================================================="
echo "Running $SDIR/postProcess.sh"
echo
echo
echo $SDIR/postProcess.sh $ASSAY $PROJNO
$SDIR/postProcess.sh $ASSAY $PROJNO

convert $PROJNO/*png $PROJNO/${PROJNO}___seqSeg.pdf
