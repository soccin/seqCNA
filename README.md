# CopyNumber analysis from sequence data using seqDNAcopy

## dev/feature/targeted-v4 branch

Copynumber analysis for sequence based assays (shallow Whole Genome sequncing, or targeted assays) using `seqDNAcopy` package from `seshanv@mskcc.org`

_N.B._ these script use a forked version of `seqDNAcopy` which has been updated to work with both human and mouse.

### getPairedCounts

```{bash}
usage: getPairedCounts.R \
	TUMOR=/path/tumor.bam \
	NORMAL=/path/normal.bam \
	ODIR=[.] \
	GENOME=[hg19] \
	GCNORM=[FALSE]
```

This script will just get the read counts for a tumor normal pair. GCNORM arg can be set to gc-normalize them.

