#!/bin/bash

#git clone git@github.com:soccin/seqCNA.git

if [ ! -e pipeline/alignments ]; then
    echo
    echo "Need to link up pipeline directory"
    echo "    ln -s /ifs/res/seq/pi/invest/Proj_No/r_00x pipeline"
    echo
    exit
fi

ls pipeline/alignments/*.bam | sed 's/.*_recal_s_/s_/' | sed 's/.bam/__/' >samples
if [ ! -e tumors ]; then
    echo
    echo "Need to define tumor and normal samples"
    echo
    exit
fi

ls pipeline/alignments/*.bam \
    | xargs -n 1 \
        bsub -o LSF.FIX/ -J FIX -n 2 -R "rusage[mem=8]" ./seqCNA/fixChromosomeNames.sh

bSync FIX

ERR1=$(parseLSF.py LSF.FIX/* | fgrep -v Succ)

if [ "$ERR1" != "" ]; then
    echo "ERROR @ Stage1"
    parseLSF.py LSF.FIX/* | fgrep -v Succ
    exit
fi

ls bamRelabel/*/*bam | fgrep -f tumors >tumorBams
ls bamRelabel/*/*bam | fgrep -f normals >normalBams

for tumor in $(cat tumorBams); do
    for normal in $(cat normalBams); do
        ./seqCNA/seqCNA.sh 100 $normal $tumor
    done
done

bSync "SEQSEG_s.*"
bSync "WGSCNA_s.*"

ERR2=$(parseLSF.py LSF/*/* | fgrep -v Succ)
if [ "$ERR1" != "" ]; then
    echo "ERROR @ Stage2"
    parseLSF.py LSF/*/*| fgrep -v Succ
    exit
fi

./seqCNA/selectBestMatch out
# Take best match from T/N Pairs
#cat pipeline/*_sample_pairing.txt | awk '{print $2"__"$1}' | fgrep -v POOL >bestMatches____out

mkdir outAll
rsync -avP --link-dest=../out out/ outAll

ls -d out/s_/*/* | fgrep -vf bestMatches____out | xargs -t rm -rf

GENOME=$(./seqCNA/GenomeData/getGenomeBuildBAM.sh $(ls pipeline/alignments/*.bam | head -1))
PROJNO=$(echo $(ls pipeline/*request.txt) | perl -ne 'm|/(Proj_.*)_request|; print $1')

if [ $GENOME == "mm10" ]; then
    ASSAY=M-IMPACT_v1
else
    echo "Unknow ASSAY" $GENOME
    exit
fi

./seqCNA/postProcess.sh $ASSAY $PROJNO
convert $PROJNO/*png $PROJNO/${PROJNO}___seqSeg.pdf