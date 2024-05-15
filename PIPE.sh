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
module load samtools
MAJOROS=7
export R_VERSION=$(R --version | head -1 | awk '{print $3}')
export R_LIBS=/home/socci/lib/R/CentOS${MAJOROS}/$R_VERSION

R --version

SDIR="$( cd "$( dirname "$0" )" && pwd )"

function usage {
    echo
    echo "   usage: $(dirname $SDIR)/PIPE.sh [sWGS|TARGETTED]"
    echo "       sWGS       Shallow Whole Genome Seq"
    echo "       TARGETTED   Targeted assays (IMPACT/EXOME/...)"
    echo
    exit
}

if [ "$#" -ne 1 ]; then
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

BAMDIR=$(ls -d pipeline/* | egrep "alignments|bam")

if [ ! -e "$BAMDIR" ]; then
    echo
    echo "Need to link up pipeline directory"
    echo "    ln -s /ifs/res/seq/pi/invest/Proj_No/r_00x pipeline"
    echo
    exit
fi

if [ ! -e samples ]; then

    ls $BAMDIR/*.bam \
        | xargs -n 1 basename \
        | sed 's/.*_recal_s_/s_/' \
        | sed 's/^Proj_.*s_/s_/' \
        | sed 's/.bam/__/' \
        >samples

fi

if [ ! -e tumors ]; then

    echo
    echo "Need to define tumor and normal samples"
    echo
    echo "And then run as bsub or nohup"
    echo
    exit

fi

echo "=============================================================================="
echo "Running $SDIR/fixChromosomeNames.sh"
echo
echo

FIXTAG=${TAG}_FIX_${UUID}
ls $BAMDIR/*.bam \
    | xargs -n 1 \
        bsub -o LSF.FIX/ -J $FIXTAG -W 59 $SDIR/fixChromosomeNames.sh

bSync $FIXTAG

echo "Done with fixChromosomeNames"
echo


ERR1=$(parseLSF.py LSF.FIX/* | fgrep -v Succ)

if [ "$ERR1" != "" ]; then
    echo "ERROR @ Stage1"
    parseLSF.py LSF.FIX/* | fgrep -v Succ
    exit
fi

ls bamRelabel/*/*bam | fgrep -f tumors >tumorBams
ls bamRelabel/*/*bam | fgrep -f normals >normalBams

nTumBams=$(wc -l tumorBams | awk '{print $1}')
nNormBams=$(wc -l normalBams | awk '{print $1}')

if [[ "$nTumBams" == "0" || "$nNormBams" == "0" ]]; then
    echo
    echo Something went wrong, No tumor/normal bams found
    echo
    exit
fi


echo "=============================================================================="
echo "Running $SDIR/seqCNA.sh"
echo
echo


for tumor in $(cat tumorBams); do
    for normal in $(cat normalBams); do
        $SDIR/seqCNA.sh $BINSIZE $normal $tumor
        EXITCODE=$?
        if [ "$EXITCODE" != "0" ]; then
            echo
            echo "ERROR in "$SDIR/seqCNA.sh $BINSIZE $normal $tumor
            echo "CODE = "$EXITCODE
            echo
            exit 1
        fi
    done
done


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
    ASSAY=M-IMPACT_v2
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
