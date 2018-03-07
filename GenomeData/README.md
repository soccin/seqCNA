# GenomeData
Tools to detect which genome was used and data files with paths / info for genomes used

Example code:
```bash
# Get Genome BUILD

GENOME_BUILD=$($SDIR/GenomeData/getGenomeBuildBAM.sh $BAM)
GENOME_SH=$SDIR/GenomeData/genomeInfo_${GENOME_BUILD}.sh
if [ ! -e "$GENOME_SH" ]; then
    echo "Unknown genome build ["${GENOME_BUILD}"]"
    exit
fi
echo "Loading genome [${GENOME_BUILD}]" $GENOME_SH
source $GENOME_SH
echo GENOME=$GENOME
```

