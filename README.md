# CopyNumber analysis from sequence data using seqDNAcopy

## Branch: master

Copynumber analysis for sequence based assays (shallow Whole Genome sequncing, or targeted assays) using `seqDNAcopy` package from `seshanv@mskcc.org`

_N.B._ these script use a forked version of `seqDNAcopy` which has been updated to work with both human and mouse.

### getPairedCounts

```
getPairedCounts [options] TUMOR=tumor.bam NORMAL=normal.bam

Compute the binned counts from a sample pair. By default normalize
for library size and GC variation.

Options:

TUMOR=File          TUMOR sample BAM file.  Required.

NORMAL=File         NORMAL sample BAM file.  Required.

GENOME=String       Genome.  Default value: hg19. Possible values: {hg19, mm10}

ODIR=Directory      Directory to write output files.  Default: "."

GCNORM=Boolean      Flag controlling whether GC-normalization is done.
                    Default value: TRUE. Possible Values: {FALSE, TRUE}

SAMPLEID=String     Sample Id. If not specified then will form an id from
                    "cleaned" version of tumor and normal bam filenames
```


### seqSegment

```
USAGE: seqSegment [options] COUNTS=counts.rda

Run the CBS algorithm on binned counts from getPairedCounts.

Options:

COUNTS=File                   Counts file created by getPairedCounts (Rdata file).  Required.

ODIR=Directory                Directory to write output files.  Default value: "."

BINSIZE=(Integer|enum)        Size of bins in base pairs. Default value: "auto". Possible
                              values: {auto,INTEGER}. If set to "auto" then determine bin
                              size dynamically
```
