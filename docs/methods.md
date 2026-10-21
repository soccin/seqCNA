# Methods

## seqCNA: copy number changes from sequencing data

seqCNA estimates copy number changes in a tumor from the ratio of its
read coverage to that of a normal sample. It uses seqDNAcopy, V.
Seshan's package of the read counting, GC correction and segmentation
steps, built on DNAcopy, with mouse support added. The values below
are the ones the pipeline uses.

### Read counts

Fragments are counted in the BAM files of the tumor and the normal. A
fragment is counted once, from the first read of the pair, if the read
passes quality control, is not marked as a duplicate, is not a
secondary alignment, has a mapped mate and a mapping quality of at
least 20, and the fragment length (insert size) is 75 to 750 bp. The
mapping quality filter is not applied on chromosome X. Each fragment
is counted at its midpoint in 100 bp bins along the genome. Fragments
whose midpoint falls in a blacklist of problematic regions are
removed (436 regions, 11.7 Mb, in the human genome; there is no
blacklist for mouse). The autosomes and chromosome X are counted;
chromosome Y and the mitochondrial genome are not.

### GC correction

The tumor counts are corrected for GC bias relative to the normal.
Each 100 bp bin is given the GC content of the 1 kb window centered on
it (pctGCdata). The tumor and normal counts are summed over the bins
of each GC value, and the log2 ratio of tumor to normal is computed
for each GC value, after scaling the tumor by the ratio of the total
normal to the total tumor counts. A loess curve (span 0.25) is fitted
to these log2 ratios against GC content, weighted by log2 of the sum
of the tumor and normal counts at each GC value. The tumor count of
each bin is then multiplied by 2^-f, where f is the fitted value at
the GC content of the bin.

### Bin size

For targeted sequencing (exome, IMPACT) the 100 bp bins are segmented
as they are. For shallow whole-genome sequencing they are merged into
larger bins, with the size chosen for each pair so that a bin at the
lower quartile of the normal counts holds about 100 normal fragments:
the normal counts are summed in 50 kb bins (500 kb if the lower
quartile is below 1), and the bin size is 100 x 50 kb divided by
their lower quartile, rounded down to a multiple of 100 bp.

### Segmentation

Bins with 35 or fewer normal fragments are removed; in targeted data
this keeps the bins in well-covered capture targets. The log2 ratio of
each remaining bin is log2(s x tumor + 1) - log2(normal + 1), where s
is the ratio of the total normal to the total tumor counts of these
bins. The log2 ratios are segmented by circular binary segmentation
(CBS, DNAcopy::segment [Olshen 2004, Venkatraman 2007]) with each bin
weighted by log2(normal + 1 - 35), the default significance level
(alpha = 0.01), and undo.splits = "sdundo" with undo.SD = 2: a change
point is removed if the means on its two sides differ by less than
two standard deviations. CBS computes its p-values by permutation; the
random seed is fixed (42), so the same counts always give the same
segments.

### Segment clustering and the diploid level

Segments at nearly the same level are merged into clusters, which
removes the small differences between segments that are at the same
copy number (adapted from FACETS [Shen 2016]). The segment means are
sorted, and the two adjacent levels that are closest are merged as
long as they differ by less than 0.04; the level of a merged cluster
is the weighted mean of the log2 ratios of its bins. Each segment then
takes the level of its cluster. The cluster with the smallest
abs(level)/sqrt(number of bins) is the diploid level and is set to 0;
all levels are shifted by the same amount, and the log2 ratios of the
bins move with their segment.

The diploid level is the level of a large part of the genome. In a
tumor with a genome doubling or many arm-level changes it may not be
two copies.

### Choice of the normal

Each tumor is compared with every normal in the normal set, which are
normal samples processed and sequenced in the same way as the tumors;
a sample that is in both sets is not compared with itself. For each
pair the root mean square of the bin log2 ratios after clustering,
sqrt(mean(log2R^2)), is computed over all bins, and the pair with the
smallest value is kept: the normal whose coverage pattern is closest
to that of the tumor. Only that pair is reported.

### Significant segments

As in [Cheng 2015], the mean and standard deviation of the bin log2
ratios in the diploid cluster are computed. For each segment a
two-sided normal test compares its level with the diploid mean, using
the diploid standard deviation as the scale. The p-values of the
segments of a tumor are adjusted for multiple testing by the
Benjamini-Hochberg false discovery rate (FDR) [Benjamini 1995]. A
segment is significant if FDR < 0.05 and the absolute log2 ratio is
greater than 1 (more than a two-fold change).

### Gene-level calls

Genes are represented by the exons of one canonical transcript per
gene. For the human exome (GRCh37/b37) these are Ensembl release 75
transcripts, with the canonical transcript of each gene from the
isoform lists of vcf2maf v1.6.12 (the MSKCC and UniProt overrides).
For the Twist mouse exome (GRCm38/mm10) they are Ensembl release 83
transcripts, with the transcript that overlaps the most capture
targets as the canonical one.

A gene is in a segment if one of its exons overlaps the segment. If a
gene overlaps more than one segment, the segment with the lowest log2
ratio is used: a deletion of part of a gene is called, and a gain of
part of a gene usually is not. A gene is called if that segment is
significant.

Calls are made on the autosomes only. A call on chromosome X would
need a sex-matched normal and the known sex of the sample, which are
rarely available. Chromosome X is still counted, segmented and
plotted.

### Software

seqCNA 6.0.0 with R 4.5.1, DNAcopy 1.84.0, seqDNAcopy 0.0.4,
pctGCdata 0.3.0, Rsamtools 2.26.0 and data.table 1.18.6.1. DNAcopy
1.75.5 fixed the weights used in weighted segmentation (earlier
releases segmented every chromosome after the first with the weights
of the first chromosome), so the results differ from those of seqCNA
5.x, which ran with an earlier DNAcopy.


## References

- Benjamini, Y, Hochberg, Y; J R Stat Soc B 1995; 57(1): 289-300.

- Cheng, DT, et al; J Mol Diagn 2015 May; 17(3): 251-264.

- Olshen, AB, et al; Biostatistics 2004 Oct; 5(4): 557-572.

- Shen, R, Seshan, VE; Nucleic Acids Res 2016; 44(16): e131.

- Venkatraman, ES, Olshen, AB; Bioinformatics 2007; 23(6): 657-663.


## Source code

- https://github.com/soccin/seqCNA
- https://github.com/soccin/seqDNAcopy
