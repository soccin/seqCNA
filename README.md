# CopyNumber analysis from sequence data using seqDNAcopy

Copy number analysis for sequence-based assays (shallow Whole Genome
sequencing, or targeted assays) using the `seqDNAcopy` package from
`seshanv@mskcc.org`.

_N.B._ these scripts use a forked version of `seqDNAcopy` which has
been updated to work with both human and mouse.

Runs on IRIS (RHEL 8, Slurm) since 6.0.0. The last JUNO/LSF version is on
branch `juno` (3e6c18d). See `CHANGELOG.md`, `docs/SLURM_PORT.md`,
`00.ISSUES.md`.

## Setup, once per checkout

```bash
module load R/4.5.1_nds
git clone --recurse-submodules git@github.com:soccin/seqCNA.git
seqCNA/00.SETUP.sh
```

`00.SETUP.sh` builds `seqDNAcopy` into `seqCNA/Rlib/`. It needs
DNAcopy, pctGCdata and RJSONIO in your R library and prints the
install commands if any is missing.

## Running

From a project directory holding `pipeline/<alignments|bam>/*.bam`
(with `.bai`), `tumors` and `normals` (one pattern per line, matched
against the BAM paths):

```bash
sbatch seqCNA/PIPE.sh TARGETTED    # IMPACT, exome: BINSIZE=100
sbatch seqCNA/PIPE.sh sWGS         # shallow WGS:   BINSIZE=auto
seqCNA/bin/checkRun.sh             # 0 ok, 1 failed, 2 running
```

Deliverables land in `<Proj_No>/`: the per-pair PNGs,
`<Proj_No>___IGV.seg`, `GeneMatrix.csv`, `SegmentMatrix.csv`,
`<Proj_No>___GeneTable.xlsx`, `<Proj_No>___seqSeg.pdf`. The best
tumor-normal pairs are in `out/`, all pairs in `outAll/`.

### Normals

- `IDT_Exome_v2_Normal`: was `~/Work/Pipelines/SeqCNV` on JUNO; not yet
  copied to IRIS (`00.ISSUES.md` #3).

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
