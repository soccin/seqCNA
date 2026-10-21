# Output

The results of a seqCNA run are delivered in the folder `seqcna/`.
`Proj_<no>` below is the project number. Each sample is named by its
tumor sample ID (`s_...`). The normal that each tumor was compared
with is chosen by the pipeline (`methods.md`, "Choice of the normal")
and is not named in the files.

| File | Content |
|---|---|
| `Proj_<no>___seqSeg.pdf` | Copy number plot of each tumor, one per page |
| `Proj_<no>___IGV.seg` | Every segment of every tumor, for IGV |
| `Proj_<no>___GeneTable.xlsx` | Called genes, one row per tumor and gene |
| `GeneMatrix.csv` | Called genes, a matrix of genes by tumors |
| `SegmentMatrix.csv` | Called genes with the segment they are in |
| `README.txt` | Instead of the three tables above, when no gene is called |
| `methods.md` | The method |
| `output.md` | This file |

## Log2 ratios

Every copy number value in the files is a log2 ratio of tumor to normal
coverage, centered so that the diploid level of the tumor is 0
(`methods.md`, "Segment clustering and the diploid level"). In a
tumor with no normal cells and a diploid genome:

| Copies | log2 ratio |
|---|---|
| 0 | far below -1 |
| 1 | -1 |
| 2 | 0 |
| 3 | 0.58 |
| 4 | 1 |
| 6 | 1.58 |
| 8 | 2 |

Normal cells in the sample move every value toward 0. With a tumor
fraction p, c copies give log2(1 - p + p x c / 2): one lost copy in a
sample with 50% tumor cells gives -0.42, not -1.

A segment is called only if its FDR is below 0.05 and its absolute
log2 ratio is greater than 1. So the calls are deletions of both
copies, in samples with more than about 50% tumor cells, and gains to
more than four copies (more than six at 50% tumor cells). Gains and
losses of one copy are not called; they can be seen in the plots and
in the IGV file.

## Proj_<no>___seqSeg.pdf

One page per tumor, titled with the tumor ID.

- The x axis is the position along the genome in Mb, with the
  chromosomes end to end. Dashed vertical lines separate the
  chromosomes, and their numbers are at the bottom. Chromosome X is
  numbered 23 (human) or 20 (mouse).
- Each point is the log2 ratio of one bin, light blue and gray on
  alternating chromosomes. The y axis runs from -5 to 5; points beyond
  that are not shown. In targeted sequencing there are bins only in the
  capture targets, so the points are not spread evenly along the
  genome.
- The red lines are the segments, at the level of their cluster.
- The dashed horizontal lines are at log2 ratios -1, 0.58 (three copies
  in a pure diploid tumor) and 1. Only segments beyond -1 or 1 can be
  called.
- The dotted horizontal lines are at plus and minus the MAD.

The two lines of text at the top are measures of the data of the
pair:

| Value | Meaning |
|---|---|
| `NumBins` | Number of bins segmented |
| `MAD` | Median absolute deviation (R `mad`, scaled to estimate a standard deviation) of the bin log2 ratios from their segment level: the noise in the data |
| `nSegs` | Number of segments, chromosome X included |
| `RMSD.Segs.noX` | Root mean square of the segment levels, chromosome X excluded; each segment counts once, whatever its length |
| `sum.logr.sq` | Sum of the squared bin log2 ratios |
| `rms.logr` | Root mean square of the bin log2 ratios. The normal with the lowest value is the one chosen |
| `rms.logr.flat` | Root mean square of the bin log2 ratios around their segment level, another measure of the noise |

Chromosome X is plotted but never called. If the tumor and the normal
come from individuals of different sex, the whole of X is near -1 (male
tumor, female normal) or near 1 (female tumor, male normal).

## Proj_<no>___IGV.seg

Every segment of every tumor, significant or not, tab-separated, one
row per segment. To view it in IGV, choose the genome of the project
(hg19 for GRCh37/b37, mm10 for mouse), then File > Load from File.
IGV colors each segment by its log2 ratio, red for a gain and blue for
a loss.

| Column | Meaning |
|---|---|
| `ID` | Tumor ID |
| `chrom` | Chromosome: 1 to 22 and X (mouse: 1 to 19 and X), without `chr` |
| `loc.start`, `loc.end` | Centers of the first and last bins of the segment, in bp |
| `num.mark` | Number of bins in the segment |
| `seg.mean` | Log2 ratio of the segment (the level of its cluster) |

The p-value and FDR of each segment are not in this file; for the
segments with a called gene they are in `SegmentMatrix.csv` and the
gene table.

## Gene call tables

The three tables list the genes called in at least one tumor
(`methods.md`, "Gene-level calls"): a gene whose segment has FDR <
0.05 and an absolute log2 ratio greater than 1, on the autosomes
only. If no gene is called in any tumor, the tables are not written,
and `README.txt` says "No significant copy number changes".

Genes are named by their gene symbol and the Ensembl ID of the
canonical transcript used.

### Proj_<no>___GeneTable.xlsx

One sheet, one row per tumor and called gene, sorted by tumor and
position.

| Column | Meaning |
|---|---|
| `ID` | Tumor ID |
| `gene`, `transcript` | Gene symbol and Ensembl transcript |
| `chrom`, `loc.start`, `loc.end` | The segment, as in the IGV file |
| `seg.mean` | Log2 ratio of the segment, to 3 decimals |
| `FDR` | FDR of the segment |

### GeneMatrix.csv

One row per called gene. The columns are `gene`, `transcript`, then
one column per tumor with at least one call. A value is the log2 ratio
of the call, or 0 if the gene is not called in that tumor. A 0 means
"not called", not "two copies": the gene can still be in a segment
that is gained or lost but not significant. Tumors without any call
have no column.

### SegmentMatrix.csv

The same calls with their segments, one row per called gene and
segment, sorted by position. A row is the call of one tumor (or of
several tumors with the same segment); its column holds the log2
ratio, and the columns of the other tumors are empty.

| Column | Meaning |
|---|---|
| `chrom`, `loc.start`, `loc.end`, `num.mark` | The segment, as in the IGV file |
| `cluster` | Number of the segment's cluster within the tumor; for checking only |
| `pValue`, `FDR` | p-value and FDR of the segment |
| `start`, `stop`, `exon` | An exon of the gene that overlaps the segment |
| `gene`, `transcript` | Gene symbol and Ensembl transcript |
| one column per tumor | Log2 ratio of the segment, in the rows of that tumor's calls |
