# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## What this is

seqCNA calls copy number changes from sequencing data (shallow WGS or
targeted/exome) by comparing tumor BAMs to normal BAMs. It uses
`seqDNAcopy`, a fork of V. Seshan's package with mouse support added.
Steps: binned read counts, GC correction, CBS segmentation, segment
clustering, FDR-based significance, gene-level calls. `METHODS.md` has
the method text used in reports; `docs/seqCNA-vs-FACETS*.md` explains
when to use this instead of FACETS.

The code is bash and R scripts for the IRIS cluster (RHEL 8, Slurm).
There is no build step, test suite, or linter. `ATTIC/` holds retired
scripts, including the LSF-era `seqCNA.sh`.

Other files to read first:
- `CHANGELOG.md`: what changed per release
- `00.ISSUES.md`: open issues
- `docs/SLURM_PORT.md`: how the port from JUNO/LSF was done and validated
- `docs/IRIS_CLUSTER_INFO.md`: cluster facts with job ids
- `00.SPEC_IRIS_PORT.md`: the port specification and decisions D1-D16

## Setup, once per checkout

```bash
git submodule update --init     # seqDNAcopy source
seqCNA/00.SETUP.sh              # builds Rlib/seqDNAcopy (gitignored)
```

R 4.5.x must be on PATH (`module load R/4.5.1_nds`). `00.SETUP.sh`
checks that DNAcopy, pctGCdata and RJSONIO load from the user's R
library and stops with install hints otherwise.

## Running

The pipeline runs from a project directory with this repo checked out
(or linked) as `seqCNA/` inside it. All outputs go to the current
directory. One run per project directory.

Required inputs in the project directory:
- `pipeline/` containing a subdirectory matching `alignments|bam` with
  the BAMs, each with a `.bai`
- `tumors`, `normals`: one pattern per line, matched with `fgrep -f`
  against BAM paths

```bash
sbatch seqCNA/PIPE.sh TARGETTED   # BINSIZE=100 (IMPACT, exome)
sbatch seqCNA/PIPE.sh sWGS        # BINSIZE=auto
seqCNA/bin/checkRun.sh            # 0 ok, 1 failed, 2 running, 3 sacct failed
```

`sbatch` is the `~/bin` wrapper; it exports `SBATCH_SCRIPT_DIR`, from
which `PIPE.sh` finds the checkout. Jobs are submitted with `~/bin/bsub`
and waited on with `~/bin/bSync.sh`.

Single steps can be run by hand. Any wrapped tool run with no arguments
prints `docs/<name>.doc`.

```bash
seqCNA/getPairedCounts GENOME=hg19 TUMOR=t.bam NORMAL=n.bam ODIR=out/x SAMPLEID=s_T__s_N
seqCNA/seqSegment COUNTS=out/x/s_T__s_N_Counts.rda ODIR=out/x BINSIZE=100
seqCNA/selectBestMatch out
seqCNA/postProcess.sh <Assay> [Proj_XXXXX]
```

Valid `<Assay>` values are the directory names under
`resources/geneAnnotations/` (excluding `src`). PIPE.sh chooses the
assay from the genome only (b37/hg19 -> Exome, mm10/GRC_m38 ->
Mouse-Twist); other assays (IMPACT410, IMPACT505, M-IMPACT_v*) require
running postProcess.sh directly.

## Pipeline flow (PIPE.sh)

PIPE.sh is the Slurm control job (`#SBATCH` header: `cmobic_cpu`,
`--qos=priority`, 3 days, log `SLURM.CTRL/<jobid>.out`). Each stage is
fanned out with `seqSub` (`bin/slurmTools.sh`) through `bin/runStage.sh`
and waited on with `waitStage` = `bSync` + `bCheck` (sacct State must
be `COMPLETED`). Job names are `qseqCNA_<STAGE>_<pid>`.

1. CHECK (control job): BAM directory, BAMs and indexes, `tumors` and
   `normals`, one genome across BAMs via `GenomeData/getGenomeBuildBAM.sh`
   and the assay from it, project id from `pipeline/*request.txt` or
   `extractProjectIDFromPath.py`. Nothing is submitted until all pass.
2. FIX, one job per BAM: `fixChromosomeNames.sh` writes
   `bamRelabel/<base>/<base>__relabel.bam` with `chr` stripped from
   contig names. b37 BAMs are symlinked with their index instead. Logs
   in `SLURM.01.FIX/`.
3. After FIX: `tumorBams`, `normalBams` from the relabeled BAMs; the
   `GENOME=` argument from the relabeled tumor's tag; the pair list
   `out/pipeline_info/pairs.tsv` (PAIR, TUMOR_BAM, NORMAL_BAM, ODIR) for
   every tumor x normal; `out/pipeline_info/cmd.sh.log`.
4. COUNT, one job per pair: `getPairedCounts`. Class from `runClass`
   on the BAM sizes. Logs in `SLURM.02.COUNT/`.
5. SEG, one job per pair: `seqSegment` on the `_Counts.rda`. Logs in
   `SLURM.03.SEG/`.
6. BEST (control job): `selectBestMatch out` keeps, for each tumor, the
   normal with the lowest `rms.logr` and writes `bestMatches____out`.
   `out/` is hard-linked to `outAll/`, then the non-best pair
   directories are deleted from `out/`.
7. POST, one job: `postProcess.sh <Assay> <ProjNo>` copies PNGs to
   `<ProjNo>/`, concatenates the `.seg` files into `<ProjNo>___IGV.seg`,
   runs `fixXChrom.R` (renames the highest numeric chromosome to X),
   `getGeneCalls` (`GeneMatrix.csv`, `SegmentMatrix.csv`,
   `<ProjNo>___GeneTable.xlsx`, or `README.txt` when nothing is
   significant), and `magick` merges the PNGs into
   `<ProjNo>___seqSeg.pdf`. Logs in `SLURM.04.POST/`.
8. DELIVERABLES and STAGING (control job): every best-pair `.seg`,
   `.out`, `.png` and every `<ProjNo>/` file must be non-empty;
   `bCheckAll` sweeps every job of the run; `STATUS=COMPLETED`.

Run records: `00.RUNSTATUS.txt` (STATUS, STAGE, RC, MESSAGE, ...),
`SLURM.CTRL/jobs.tsv` (job id, stage, sample, log), `#SEQCNA_EXIT=<rc>`
as the last line of every stage log. Any exit before COMPLETED marks
the run FAILED and cancels its queued jobs (EXIT trap).

After the pipeline:
- `delivery.sh /path/Proj_X/r_00N` copies deliverables into `seqcna/`
  there (not ported; JUNO paths). `DELIVERY.md` and `DELIVERY_ZIP.md`
  are the email templates.
- `Rscript --no-init-file seqCNA/redactNormals.R` removes `__<normal>`
  from sample IDs in the four deliverable tables, rewriting them in
  place. Run it from a directory that contains exactly one copy of each.
- `Rscript --no-init-file seqCNA/plotSeg.R [-r] SEG.rda` re-plots a
  `_seqSeg.rda`; `-r` removes the normal ID from the title.

## Architecture

### Wrapper pattern

`getPairedCounts`, `seqSegment`, `getGeneCalls`, and `selectBestMatch`
are symlinks to `execWrapper.sh`. The wrapper exports `SDIR` (repo
directory) and `SNAME` (invoked name), requires R 4.5.x on PATH and
`Rlib/seqDNAcopy`, prepends `$SDIR/Rlib` to `R_LIBS`, and runs
`Rscript --vanilla $0.R "$@"`. To add a tool: write `foo.R`, symlink
`foo -> execWrapper.sh`, add `docs/foo.doc`.

- `--vanilla` skips `~/.Rprofile` and `~/.Renviron` but still honors
  the `R_LIBS` and `R_LIBS_USER` environment variables, which is how
  the user's library is found. Wrapped scripts cannot use `~/.Rprofile`
  helpers; `include/tools.R` redefines `cc()` and `len()`.
- Wrapped scripts source `include/*.R` via `Sys.getenv("SDIR")`, which
  falls back to the relative path `seqCNA`.
- Standalone scripts (`plotSeg.R`, `redactNormals.R`, `fixXChrom.R`)
  run with `Rscript --no-init-file` and must not depend on `~/.Rprofile`.

### Shell helpers (bin/)

- `loadTools.sh`: sourced by `PIPE.sh`, `fixChromosomeNames.sh`,
  `postProcess.sh`. Loads `samtools/1.20`, prepends ImageMagick 7,
  checks R 4.5.x, sets `TMPDIR` under `/localscratch/$USER`.
- `slurmTools.sh`: `bSync`, `bCheck`, `bCheckAll`, `bKill`,
  `failedJobsMessage`. ChIP-seq's file with `CHIP_` -> `SEQCNA_`; keep
  them in step.
- `runStage.sh`: per-job `TMPDIR`, exit trailer. `checkRun.sh`: verdict
  from the status file, sacct and the trailers.

### Argument parsing

Wrapped scripts take `KEY=VALUE` arguments parsed by
`include/parseArgs.R`. Keys are uppercased. Values are coerced to the
class of the default. A `NULL` default marks a required argument; if
one is missing, `parseArgs` returns NULL and the script prints its doc
file.

### R and seqDNAcopy

- One R, 4.5.x, for every script. `seqDNAcopy/` is a git submodule with
  the fork's source (`git submodule update --init`). `00.SETUP.sh`
  builds it into `Rlib/` (gitignored) and records versions in
  `Rlib/VERSIONS.txt`. The fork's NAMESPACE imports `IRangesList`
  (0.0.4); `RangesList` is gone from IRanges 2.x.
- DNAcopy, pctGCdata and RJSONIO are not in the shared R install; they
  live in the user's library.

### Genome handling

`GenomeData/getGenomeBuildBAM.sh` identifies the reference by the md5
of the sorted `@SQ` header lines (first 3 columns); it needs `samtools`
on PATH and exits 1 if the header cannot be read. The resulting tag is
mapped in `PIPE.sh` in two places that must stay consistent with it:
- before FIX, on the original BAM: tag -> annotation assay (Exome or
  Mouse-Twist)
- after FIX, on the relabeled tumor BAM: tag -> `hg19` or `mm10` (the
  `GENOME=` argument). The relabeled md5 differs from the original's
  (UCSC mm10 becomes `mm10-relabel`).
`getGeneCalls.R` maps hg19 to b37 (`normalizeGenomeTag`) when choosing
`resources/geneAnnotations/<ASSAY>/<genome>/gene_annotations.txt.gz`.
Chromosome X is numeric in segment output (23 human, 20 mouse).

### Sample ID conventions

`getSampleTag.sh` reads `SM:` from the BAM `@RG` lines and adds an `s_`
prefix if missing. A pair ID is `<tumor>__<normal>`. Several scripts
parse this format with regexes: `selectBestMatch.R` (`__s_`),
`getGeneCalls.R` (sample columns match `^s_`), `redactNormals.R`
(`__.*`), `plotSeg.R` (`__s_.*$`). Changing the naming breaks them.

### Per-pair outputs

`out/<first 2 chars of tumor ID>/<tumorId>/<pairId>/` contains:
- `<pairId>_Counts.rda` and `.json` (getPairedCounts)
- `<pairId>_seqSeg.rda` (DNAcopy object named `out`), `.seg`, `.png`,
  `.out` (seqSegment)

`seqsegment()` returns positions in Mb; seqSegment.R multiplies by 1e6
before writing `.seg`. The `.out` file is `key = value` QC metrics.
`parseSeqCNAOutFile()` in `include/tools.R` reads it and treats a file
without `X.seg.mean.avg` as incomplete. selectBestMatch uses it to rank
normals; getGeneCalls uses it to get the genome. Outputs record
`VERSION` and the repo's `git describe` (suffix `-UNCOMMITED` if the
tree is dirty).

### Calling parameters

- seqSegment.R: `MINBINCOUNT=35`, `undo.SD=2`, cluster merge threshold
  0.04 (`include/clusterSegs.R`, adapted from FACETS). The cluster with
  the smallest `abs(logR)/sqrt(nProbes)` is set to logR 0.
- getGeneCalls.R: significant if `FDR < 0.05` and `|seg.mean| > 1`.
  When a gene overlaps several segments, the lowest `seg.mean` is used,
  so deletions take priority over amplifications. For the Exome assay
  the canonical isoforms come from `resources/vcf2maf/v1.6.12/`.
- `METHODS.md` states these values; update it if they change.

## Environment dependencies

- `~/bin/sbatch`, `~/bin/bsub`, `~/bin/bSync.sh`: the personal Slurm
  wrappers, as in ATAC-seq and ChIP-seq. `docs/IRIS_CLUSTER_INFO.md`.
- `module load R/4.5.1_nds` (personal module in `~/.modules`) before
  submitting; batch jobs inherit the environment (`--export=ALL`).
- `samtools/1.20` module; ImageMagick at
  `/usersoftware/core001/common/RHEL_8/imagemagick/7.1.2-33`.
- `/localscratch/$USER` on every node for temp files.
