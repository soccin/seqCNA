# sWGS CopyNumber analysis using seqDNAcopy

Copynumber analysis for shallow Whole Genome sequncing data using `seqDNAcopy` package from `seshanv@mskcc.org`

_N.B._ this script uses a forked version of `seqDNAcopy` which has been updated to work with both human and mouse.

```{bash}
usage: doWGSCNA.R \
	TUMOR=/path/tumor.bam \
	NORMAL=/path/normal.bam \
	BINSIZE=[auto] \
	MAPLOC=[false] \
	GENOME=[hg19]
```

