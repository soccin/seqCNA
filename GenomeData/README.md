# GenomeData

`getGenomeBuildBAM.sh BAM` prints a genome tag for a BAM from the md5 of
its sorted `@SQ` header lines (first three columns), or `unknown <md5>`.
It needs `samtools` on PATH. `getGenomeDictMD5.sh` computes the same md5
from a header on stdin, for adding a new genome to the `case` table.

```bash
GENOME_BUILD=$($SDIR/GenomeData/getGenomeBuildBAM.sh $BAM)
```

The tag is mapped to the `GENOME=` argument of `getPairedCounts`
(hg19 or mm10) and to the gene annotation assay in `PIPE.sh`. Keep the
three places consistent. The old `genomeInfo_<tag>.sh` files, which held
`/ifs` paths and were never sourced by the pipeline, are in
`attic/GenomeData/`.
