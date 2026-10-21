# seqCNA Pipeline: Description and Comparison with FACETS

## Overview

The seqCNA pipeline is an in-house bioinformatics tool for somatic copy number analysis
from tumor/normal sequencing data. It implements a "total copy number" approach, in contrast
to allele-specific methods such as FACETS. The pipeline was developed at the MSK
Bioinformatics Core and is built around `seqDNAcopy`, V. Seshan's package for read
counting, GC correction and segmentation. The full method is in `docs/methods.md`.

Source code: https://github.com/soccin/seqCNA

---

## What seqCNA Does

### Algorithm

1. **Read counting and binning.** For each tumor/normal pair, the pipeline counts
   fragments in 100 bp bins, each fragment once at its midpoint. For targeted/exome data
   (MSK-IMPACT, WES) these bins are segmented as they are; for shallow WGS they are merged
   into larger bins, with the size chosen for each pair from the normal's coverage.

2. **GC normalization.** Raw bin counts are corrected for GC-content bias, which is the
   dominant technical confounder in coverage-based copy number analysis.

3. **Log-ratio computation.** For each bin the log2 ratio of tumor-to-normal coverage is
   computed. This is the central signal.

4. **Segmentation via CBS.** The log-ratio track is segmented using the Circular Binary
   Segmentation (CBS) algorithm [Olshen 2004] as implemented in `seqDNAcopy`. CBS
   identifies statistically significant breakpoints where copy number changes.

5. **Clustering and baseline.** Segments at nearly the same level are merged into
   clusters. The cluster with a level near 0 and many bins is the baseline and is set to
   0, so positive values are gains and negative values losses relative to the tumor's
   baseline level. Each segment gets a p-value against the bins of the baseline cluster,
   following the significance test of the MSK-IMPACT method [Cheng 2015], and an FDR
   over the segments of the tumor.

6. **Gene-level calls.** Each gene is given the log2 ratio of one segment: of the
   segments its exons overlap, the one with the lowest log2 ratio. The gene is called if
   that segment has FDR < 0.05 and |log2R| > 1, on the autosomes only. A call is a log2
   ratio; it is not put into a category such as gain, amplification or deletion.

### Inputs

- Tumor BAM files: hg19/b37 or mm10 (mouse support via a forked version of `seqDNAcopy`)
- Normal BAM files: matched, unmatched or pooled. Each tumor is compared with every
  normal, and the normal whose coverage pattern is closest to the tumor's (lowest
  `rms.logr`) is kept

### Outputs

Described in `docs/output.md`, which is delivered with the results.

| File | Description |
|------|-------------|
| `Proj_<no>___seqSeg.pdf` | Copy number plot of each tumor, one per page |
| `Proj_<no>___IGV.seg` | Every segment of every tumor, for IGV |
| `Proj_<no>___GeneTable.xlsx` | Called genes, one row per tumor and gene, with the FDR |
| `GeneMatrix.csv` | Called genes, a matrix of genes by tumors |
| `SegmentMatrix.csv` | Called genes with their segments, p-values and FDR |

---

## What FACETS Does

FACETS (Fraction and Allele-specific Copy number Estimate from Tumor/Normal Sequencing)
[Shen & Seshan 2022] is a more sophisticated algorithm that models allele-specific copy
number jointly with tumor purity and ploidy. It performs joint segmentation of the total
copy number log-ratio and the allele-specific signal (B-allele frequency from germline
heterozygous SNPs). From this it infers:

- Integer total copy number (TCN) per segment
- Minor (lesser) copy number (LCN) per segment, enabling LOH calls
- Cellular fraction (tumor purity)
- Ploidy estimate
- Copy-neutral LOH, heterozygous deletions, homozygous deletions, allele-specific gains

FACETS is the standard copy number method in TEMPO (the MSK clinical/CMO production
pipeline) and is widely used for published cancer genomics studies.

---

## Comparison: seqCNA vs FACETS

### Conceptual differences

| Property | seqCNA | FACETS |
|----------|--------|--------|
| Copy number model | Total (log-ratio only) | Allele-specific (integer TCN + LCN) |
| Tumor purity/ploidy | Not modeled | Explicitly estimated |
| LOH detection | No | Yes |
| Output units | Continuous log2 ratio | Integer copy number |
| Shallow vs deep deletion | Cannot distinguish | Distinguishes via LCN=0 |
| Intratumor heterogeneity | Not modeled | Partially modeled (cf per segment) |

### Practical/operational differences

| Property | seqCNA | FACETS |
|----------|--------|--------|
| Normal requirement | Matched, unmatched, or pooled | Must be matched tumor/normal |
| Mixed bait sets | Tolerant | Cannot run (different capture arrays incompatible) |
| Sample types | sWGS, IMPACT, WES, mouse IMPACT | WGS, WES, targeted panels |
| Mouse genome support | Yes (mm10 via forked seqDNAcopy) | No (human only) |
| Parameter complexity | Low (few key parameters) | Higher (cval, dipLogR, purity init, etc.) |
| Failure rate | Low; degrades gracefully | Non-trivial failure rate on low-purity or poor-QC samples |
| Segmentation on WGS | Appropriate | Prone to hypersegmentation (bin size needs tuning) |

### When to use seqCNA

- Samples lack a matched normal (unmatched or pooled-normal run required)
- Shallow whole-genome sequencing data
- Mouse tumor samples
- Cross-platform cohorts (e.g., some samples on IMPACT, some on WGS) where a single
  consistent method is preferable over running two different algorithms
- Results that follow the significance test of the clinical MSK-IMPACT copy number
  method [Cheng 2015]

### When to use FACETS

- Matched tumor/normal pairs with sufficient sequencing depth
- Allele-specific copy number needed (e.g., distinguishing LOH from homozygous deletion)
- Integer copy number calls required (e.g., for clinical reporting, OncoPrint with
  Shallow/Deep/Gain/Amplification categories)
- Tumor purity and ploidy estimates are scientifically important
- Integration with TEMPO or downstream tools that consume FACETS output format

### Known limitations of seqCNA

**Fundamental methodological constraint: relative copy number only**

seqCNA is a relative copy number method. Its primary output is the tumor/normal log2
ratio per bin/segment. It does not model, estimate, or output:

- Absolute (integer) copy number
- Tumor purity (cellular fraction)
- Ploidy
- Loss of heterozygosity (LOH)
- Allele-specific copy number

This is not an implementation gap but a design-level constraint: the pipeline uses only
read-depth and does not incorporate the B-allele frequency of germline heterozygous
SNPs. The log2 ratios are relative to a baseline level of the tumor, which is not known
to be two copies. Without purity and ploidy a log2 ratio cannot be converted to a number
of copies, and the method provides no way to estimate them. For studies where
absolute copy number, purity estimates, or LOH calls are scientifically required, seqCNA
is not a substitute for FACETS; it is explicitly a fallback for situations where FACETS
cannot be run.

**Call thresholds are fixed, not data-driven**

Significant gene-level calls are reported when FDR < 0.05 and |log2R| > 1 (i.e., a
coverage ratio more than twice or less than half that of the tumor's baseline cluster).
This threshold is hard-coded and not adapted to tumor purity. As purity drops, the log2
ratio of every change moves toward 0, so real events can fall under the fixed threshold
and be missed. Conversely, noisy samples
may generate many false positives near the threshold. There is no purity-corrected
significance model as there is in FACETS.

**Cannot distinguish shallow from deep deletion**

Because only total depth is available and purity is not estimated, the depth of a
deletion cannot be read from its log2 ratio: the value does not say whether one copy or
both copies are lost. The user has to make a binary choice about how to label deletions
in downstream tools (e.g., oncoprint). This is less of a problem for deep
focal deletions of known tumor suppressors (CDKN2A, PTEN, etc.) where the log-ratio
signal is extreme, but is a genuine ambiguity for broad, shallow losses.

**Performance on deep WGS is suboptimal**

The pipeline was developed and tuned for shallow WGS and targeted panels (MSK-IMPACT).
When applied to deep whole-genome sequencing data for the first time, both
hypersegmentation and a poorly chosen baseline cluster were observed across a subset of
samples, to the extent that gene-level calls could not be trusted without manual QC and
parameter adjustment. The high read depth at deep WGS coverage provides strong
statistical power to detect very small copy number differences that may reflect technical
noise or normal copy number variation rather than somatic events, driving oversegmentation.

**Pooled-normal caveats**

One of seqCNA's practical advantages --the ability to run without a matched normal --also
introduces a source of systematic error. When a pooled or unmatched normal is used:

- Any germline copy number variants (CNVs) in the tumor sample that are absent from
  the pooled normal will appear as somatic events. Conversely, germline CNVs in the
  pool not present in the tumor will produce spurious apparent losses.
- The normal samples in the pool must have been prepared and sequenced under the same
  library and assay conditions as the tumor. Depth differences, capture kit version
  differences, or batch effects between the tumor and the pooled normal will introduce
  systematic biases that inflate the apparent copy number signal.
- For targeted panels, if a gene is present in the tumor's capture kit but absent from
  the pooled normal's kit (or vice versa), that gene is uncallable; the coverage ratio
  is meaningless.

The pooled normal does not correct for germline genetics --only for technical/library
biases. It is appropriate when a matched normal is unavailable and a tumor-only analysis
is the only option, but results should be interpreted with extra caution.

**Reduced call rate compared to FACETS on targeted data**

When applied to the same IMPACT or exome data, seqCNA tends to produce fewer copy number
calls than FACETS. While this partly reflects the fact that seqCNA does not model purity
(so attenuated events fall below the fixed threshold), it likely also reflects genuinely
lower sensitivity for focal or low-amplitude events. Users comparing seqCNA results to
FACETS results from the same cohort should expect systematic differences in call counts.

**Intended role**

seqCNA should be understood as the method of choice for shallow WGS and for situations
where FACETS cannot run (no matched normal, mouse samples, cross-platform cohorts). For
matched human tumor/normal data from IMPACT or exome sequencing where sufficient depth
and heterozygous SNP density exist, FACETS is the more appropriate primary method.

### Known limitations of FACETS

**Normal sample requirements**

- Requires a matched tumor/normal pair; cannot be run on unmatched samples. For mouse
  models, the normal must come from the same individual (e.g., a tail biopsy) --a
  pooled or unrelated normal is generally insufficient.
- Cannot be run when tumor and normal were captured with different bait sets (e.g.,
  Agilent SureSelect vs. IDT v2 exome arrays). Genes present in one kit but not the
  other become uncallable, and the coverage profiles are not comparable.

**Purity and diploid-state estimation failures**

FACETS derives copy number by jointly fitting the total log-ratio signal and the
B-allele frequency (BAF) of germline heterozygous SNPs. This makes the quality of the
purity/ploidy solution highly sensitive to data quality:

- *Low tumor purity.* When the tumor cellular fraction is low, the allele-frequency
  shift at heterozygous SNPs is small and the copy number log-ratio signal is
  attenuated. FACETS may fail to identify a clean purity/ploidy solution, or may
  converge on an incorrect one.

- *Incorrect diploid baseline (dipLogR).* FACETS must estimate the log-ratio value
  corresponding to the diploid state from the data itself. In highly aneuploid tumors
  --particularly those with whole-genome duplication --there may be no large diploid
  region to anchor this estimate, causing the baseline to be set incorrectly. All
  downstream integer copy number calls then inherit this error. The FACETS authors
  acknowledge explicitly that their flagging heuristics do not cover all scenarios
  where dipLogR is unreliable: "I am not claiming that I have covered all possible
  scenarios that lead to bad dipLogR estimate." In practice, the output should always
  be manually inspected; automated flags catch only the most obvious failures.

- *Specific dipLogR warning flags.* FACETS emits two named warnings when it suspects
  a bad diploid estimate: `"mafR not sufficiently small"` (too few segments with
  balanced alleles, undermining the baseline fit) and `"could be polyclonal 1 copy
  loss"` (the data are consistent with genome doubling, meaning a diploid-looking
  segment may actually be tetraploid).

- *Highly aneuploid or genome-doubled tumors.* Whole-genome duplication is common in
  carcinomas. In a doubled genome the majority of the genome is at 4 copies, making
  it difficult for FACETS to correctly anchor the diploid state without additional
  prior information.

- *Inbred mouse models.* FACETS requires a sufficient density of germline heterozygous
  SNPs to compute BAF-based segmentation. Inbred mouse strains are nearly isogenic,
  providing too few heterozygous sites for the algorithm to work. This is a fundamental
  constraint, not a tuning issue.

**Parameter sensitivity and tuning burden**

- The `cval` parameter (segmentation penalty, analogous to a p-value threshold)
  directly controls the number and size of detected segments. There is no universally
  correct value: lower values oversegment, higher values undersegment. The appropriate
  setting differs between IMPACT, exome, and WGS data and sometimes between individual
  samples. Expert input (from the algorithm's author or experienced users) is often
  needed to review and re-tune problematic cases.

- On deep WGS data, the default FACETS pipeline (tuned for IMPACT/exome) is prone to
  hypersegmentation --generating an implausibly large number of short segments. Bin
  size and `cval` both need adjustment, and even then results may be noisier than on
  exome data.

**Operational limitations**

- Fingerprinting mismatches between tumor and normal samples will cause FACETS to
  produce nonsensical results (the allele-frequency signal reflects the germline
  difference rather than somatic change). These failures may not be immediately
  obvious from the plots without careful review.
- No native support for mouse (non-human) genomes.
- The clonal subclone analysis function (`emcncf2`) is noted in the package
  documentation as being actively reworked and should not be used for routine analysis.
- Multi-patient cohort output organization is less convenient; PDFs are not
  automatically compiled per-project as in seqCNA.

**Role of seqCNA as a fallback**

In practice, seqCNA serves as the fallback for cases where FACETS fails or cannot be
run. Its primary advantage in this role is that it follows the significance test of the
copy number method of the MSK clinical assay (DMP/MSK-IMPACT) [Cheng 2015], so its
output can be compared with clinical IMPACT results even when the more sophisticated
FACETS solution is unavailable. It is not an exact copy of that method: for example,
seqCNA merges segments into clusters with a threshold of 0.04, where DMP uses 0.08.

---

## Methods Text (for publications)

The methods text is `docs/methods.md`; it is delivered with the results of every run.

## References

- Cheng, DT, et al; J Mol Diagn 2015 May; 17(3): 251-264.

- Olshen, AB, et al; Biostatistics 2004 Oct; 5(4): 557-572.

