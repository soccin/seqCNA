# CopyNumber analysis from sequence data using seqDNAcopy

Copy number analysis for sequence-based assays (shallow Whole Genome
sequencing, or targeted assays) using the `seqDNAcopy` package from
`seshanv@mskcc.org`.

_N.B._ these scripts use a forked version of `seqDNAcopy` which has
been updated to work with both human and mouse.

Runs on IRIS (RHEL 8, Slurm) since 6.0.0. The last JUNO/LSF version is on
branch `juno` and tag `v5-juno` (3e6c18d). See `CHANGELOG.md`,
`docs/SLURM_PORT.md`, `00.ISSUES.md`.

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
(with `.bai`), first run it directly, without `tumors` or `normals`:

```bash
seqCNA/PIPE.sh TARGETTED           # writes samples and stops
```

`samples` has one line per BAM (`s_<id>__`). Put the tumor lines in
`tumors` and the normal lines in `normals` (one pattern per line,
matched against the relabeled BAM paths
`bamRelabel/<base>/<base>__relabel.bam`), then submit:

```bash
sbatch seqCNA/PIPE.sh TARGETTED    # IMPACT, exome: BINSIZE=100
sbatch seqCNA/PIPE.sh sWGS         # shallow WGS:   BINSIZE=auto
seqCNA/bin/checkRun.sh             # 0 ok, 1 failed, 2 running
```

Or run it directly, without `sbatch`, in a shell that stays up for the
whole run (an interactive job, tmux, nohup):
`seqCNA/PIPE.sh TARGETTED`. Run `checkRun.sh` on the same host.

Deliverables land in `<Proj_No>/`: `<Proj_No>___IGV.seg`,
`GeneMatrix.csv`, `SegmentMatrix.csv`, `<Proj_No>___GeneTable.xlsx`
(or `README.txt` when no gene is significant) and
`<Proj_No>___seqSeg.pdf`. They name each sample by its tumor ID only,
so they do not show which normal was used. `<Proj_No>/` also holds the
per-pair PNGs, titled with the pair ID `<tumor>__<normal>`; they are
for checking, not for delivery. The normal chosen for each tumor is in
`pairing.tsv` in the run directory. The best tumor-normal pairs are in
`out/`, all pairs in `outAll/`.

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

MINBINCOUNT=(Integer)         Minimum number of counts per bin to be used in segmentation.
                              Default value: 35.

YLIM=(Number)                 The plots show log2 ratios from -YLIM to YLIM.
                              Default value: 5.

TUMORID=String                If set, also write <sampleId>_seqSeg.redacted.png: the same
                              plot titled with this id instead of the pair id, so it does
                              not show the normal. PIPE.sh sets it to the tumor id; the
                              deliverable PDF is built from these plots. Default value: "".
```
