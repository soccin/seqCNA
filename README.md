# CopyNumber analysis from sequence data using seqDNAcopy

# dev/feature/targeted-v4 branch

Copynumber analysis for sequence based assays (shallow Whole Genome sequncing, or targeted assays) using `seqDNAcopy` package from `seshanv@mskcc.org`

_N.B._ this script uses a forked version of `seqDNAcopy` which has been updated to work with both human and mouse.

```{bash}
usage: seqCNA.R \
	TUMOR=/path/tumor.bam \
	NORMAL=/path/normal.bam \
	BINSIZE=[auto] \
	MAPLOC=[false] \
	GENOME=[hg19]
```

