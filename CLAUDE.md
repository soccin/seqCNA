# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## What this is

seqCNA calls copy number changes from sequencing data (shallow WGS or
targeted/exome) by comparing tumor BAMs to normal BAMs. It uses
`seqDNAcopy`, a fork of V. Seshan's package with mouse support added.
Steps: binned read counts, GC correction, CBS segmentation, segment
clustering, FDR-based significance, gene-level calls.

`docs/` holds the user-facing docs only: `methods.md` (the method text
used in reports) and `output.md` (each deliverable, its columns, how to
read the plots), both delivered with the results;
`seqCNA-vs-FACETS*.md` (when to use this instead of FACETS) and
`usage/<tool>.doc` (the usage text of the wrapped tools). Docs for
developers go in `docs/dev/`.

The code is bash and R scripts for the IRIS cluster (RHEL 8, Slurm).
There is no build step, test suite, or linter. `attic/` holds retired
scripts, including the LSF-era `seqCNA.sh`.

Other files to read first:
- `CHANGELOG.md`: what changed per release
- `00.ISSUES.md`: open issues
- `docs/dev/SLURM_PORT.md`: how the port from JUNO/LSF was done and
  validated
- `docs/dev/IRIS_CLUSTER_INFO.md`: cluster facts with job ids
- `attic/docs/00.SPEC_IRIS_PORT.md`: the port specification and
  decisions D1-D16

## Setup, once per checkout

```bash
git submodule update --init     # seqDNAcopy source
seqCNA/00.SETUP.sh              # builds seqDNAcopy into Rlib/ (gitignored)
```

R 4.5.1 must be on PATH (`module load R/4.5.1_nds`). `00.SETUP.sh`
builds seqDNAcopy from the submodule into `Rlib/`; every other package
comes from the shared or personal library. The packages that set the
output must be at the versions in `Rlib.lock.tsv` (see "Pinned
versions").

## Running

The pipeline runs from a project directory with this repo checked out
(or linked) as `seqCNA/` inside it. All outputs go to the current
directory. One run per project directory, except that a run that
failed in CHECK (nothing submitted) can be fixed and submitted again
there.

Required inputs in the project directory:
- `pipeline/` containing a subdirectory matching `alignments|bam` with
  the BAMs, each with a `.bai`
- `tumors`, `normals`: one pattern per line, matched with `fgrep -f`
  against the relabeled BAM paths `bamRelabel/<base>/<base>__relabel.bam`

Without `tumors` or `normals`, `PIPE.sh` writes `samples` (one line per
BAM, the base name minus any `Proj_..._s_` prefix plus `__`) and exits
0 before any run record is written. That is the first step of a new
project: run `seqCNA/PIPE.sh TARGETTED` directly, copy the lines of
`samples` into `tumors` and `normals`, then submit. The trailing `__`
matches the `__relabel` suffix, so `s_X__` does not select `s_X1`.

```bash
sbatch seqCNA/PIPE.sh TARGETTED   # BINSIZE=100 (IMPACT, exome)
sbatch seqCNA/PIPE.sh sWGS        # BINSIZE=auto
seqCNA/bin/checkRun.sh            # 0 ok, 1 failed, 2 running, 3 sacct failed
```

`sbatch` is the `~/bin` wrapper; it exports `SBATCH_SCRIPT_DIR`, from
which `PIPE.sh` finds the checkout. Jobs are submitted with `~/bin/bsub`
and waited on with `~/bin/bSync.sh`.

The full run can also be started directly (`seqCNA/PIPE.sh TARGETTED`)
in a shell that outlives it: an interactive job, tmux, nohup. `PIPE.sh`
sets `RUN_MODE=sbatch` only when `$0` is Slurm's spool copy
`job$SLURM_JOB_ID/slurm_script`; otherwise `RUN_MODE=direct`, the
process itself is the control job, `CTRL_JOBID=none`, and `checkRun.sh`
checks that `PID` on `HOST` is still `PIPE.sh` (run it on that host).
`SLURM_JOB_ID` alone cannot tell the two apart: inside an interactive
job it names the session.

Single steps can be run by hand. `getPairedCounts`, `seqSegment` and
`getGeneCalls` run with no arguments print `docs/usage/<name>.doc`;
`selectBestMatch` has no doc file and prints a one-line usage.

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

0. First step, before any run record: write `samples` if absent; stop
   with exit 0 if `tumors` or `normals` is missing or empty.
1. CHECK (control job): BAM directory, BAMs and indexes, `tumors` and
   `normals` each match a relabeled BAM path, one genome across BAMs via `GenomeData/getGenomeBuildBAM.sh`
   and the assay from it, project id from `pipeline/*request.txt` or
   `extractProjectIDFromPath.py`. Nothing is submitted until all pass.
2. FIX, one job per BAM: `fixChromosomeNames.sh` writes
   `bamRelabel/<base>/<base>__relabel.bam` with `chr` stripped from
   contig names. b37 BAMs are symlinked with their index instead. Logs
   in `SLURM.01.FIX/`.
3. After FIX: `tumorBams`, `normalBams` from the relabeled BAMs; the
   `GENOME=` argument from the relabeled tumor's tag; the pair list
   `out/pipeline_info/pairs.tsv` (PAIR, TUMOR_BAM, NORMAL_BAM, ODIR,
   TUMOR, NORMAL) for every tumor x normal;
   `out/pipeline_info/cmd.sh.log`. A sample in both `tumors` and
   `normals` is allowed; its self-pair is skipped, and every tumor needs
   one normal other than itself.
4. COUNT, one job per pair: `getPairedCounts`. Class from `runClass`
   on the BAM sizes. Logs in `SLURM.02.COUNT/`.
5. SEG, one job per pair: `seqSegment` on the `_Counts.rda`, with
   `TUMORID=<tumor>` so that it also writes the plot titled with the
   tumor ID only (`_seqSeg.redacted.png`). Logs in `SLURM.03.SEG/`.
6. BEST (control job): `selectBestMatch out` keeps, for each tumor in
   `pairs.tsv`, the normal with the lowest `rms.logr` (one per tumor; a
   tie goes to the first pair ID) and writes `bestMatches____out` and
   `pairing.tsv` (TUMOR, NORMAL, PAIR, rms.logr). `pairing.tsv` names
   the normals: it stays in the run directory and is not delivered.
   `out/` is hard-linked to `outAll/`, then the non-best pair
   directories are deleted from `out/` (whole-line match against
   `bestMatches____out`).
7. POST, one job: `postProcess.sh <Assay> <ProjNo>` copies the PNGs
   (pair IDs, for checking; not delivered) to `<ProjNo>/`, concatenates
   the `.seg` files into `<ProjNo>___IGV.seg`, runs `fixXChrom.R`
   (renames the highest numeric chromosome to X), `getGeneCalls`
   (`GeneMatrix.csv`, `SegmentMatrix.csv`, `<ProjNo>___GeneTable.xlsx`,
   or `README.txt` when nothing is significant), then
   `redactDeliverables.R`, which replaces each pair ID by its tumor ID
   in the IGV `.seg`, the two CSVs and the xlsx, using `pairing.tsv`.
   `magick` merges the `_seqSeg.redacted.png` plots into
   `<ProjNo>___seqSeg.pdf`. So every deliverable except the PNGs names
   a sample by its tumor ID only. Logs in `SLURM.04.POST/`.
8. DELIVERABLES and STAGING (control job): every best-pair `.seg`,
   `.out`, `.png`, `.redacted.png`, `pairing.tsv` and every `<ProjNo>/`
   file must be non-empty; `bCheckAll` sweeps every job of the run;
   `STATUS=COMPLETED`.

Run records: `00.RUNSTATUS.txt` (STATUS, STAGE, RC, MESSAGE,
RUN_MODE, CTRL_JOBID, HOST, PID, ...), `SLURM.CTRL/jobs.tsv` (job id,
stage, sample, log), `#SEQCNA_EXIT=<rc>` as the last line of every
stage log. Any exit before COMPLETED marks the run FAILED and cancels
its queued jobs (EXIT trap; also on SIGTERM, SIGINT, SIGHUP).

The client-facing version is static: the `- **Version:** vX.Y.Z`
line of `VERSION.md`, changed by hand only when a release is made,
with a git tag of the same name (`v6.0.0`). Never derive it from git.
`PIPE.sh` reads it with `bin/getRelease.sh`, records it as `RELEASE=`
in `00.RUNSTATUS.txt` and `cmd.sh.log`, and stops if it is missing;
`VERSION=` there is `git describe --tags --always --long --dirty`, the
exact code. `delivery.sh` puts the release of its own checkout in the
email and the `docs/output.md` and `docs/methods.md` links.

After the pipeline:
- `seqCNA/delivery.sh [-n] -d|DELIVERY_DIR[/r_NNN]`, run from the
  analysis directory, delivers a run whose `bin/checkRun.sh` exits 0
  (README.md, "Delivery"). It is ATAC-seq's `deliverResults.sh`: every
  write to the delivery area goes through `asbot` (`dzdo -u
  svc_core001_bot01`), `-d` derives the folder
  `/data1/core002/res/bic/results/<pi>/<invest>/Proj_<no>` from the
  `.../Users/<Pi>/<Invest>/Proj_<no>` path, the run folder is the next
  `r_NNN`, and `-n` prints those writes and the permissions-database
  call instead of running them. It copies the IGV `.seg`, the PDF and
  the three call tables (or `README.txt`), plus `docs/methods.md` and
  `docs/output.md`, to `r_NNN/seqcna/`; never the PNGs or `pairing.tsv`. Before copying it
  stops if a table still contains `__<normal>` for a normal in
  `pairing.tsv` (the xlsx is read with `unzip -p`). `PROJNO` and
  `GENOME_TAG` come from `00.RUNSTATUS.txt`; the tag is mapped to the
  assemblies the database accepts (b37, hg19, hybrid, mm10).
  `bin/readme2yaml.R` writes `project.yaml` (pipeline `seqcna`) from
  the first `README.txt` going up to `Proj_<no>`; there is no fallback
  to `README.md`. The email is built from `tmpldeliveryEmail_01.txt`,
  then `_02` (call tables) or `_03` (no significant calls), then `_99`,
  with `{PROJNO}`, `{RUNDIR}` and `{VERSION}` filled in.
  `DELIVERY_ZIP.md` is the email to use, by hand, when the results go
  out as a ZIP attachment.
- `seqCNA/bin/rscript.sh seqCNA/plotSeg.R [-r] SEG.rda` re-plots a
  `_seqSeg.rda`; `-r` removes the normal ID from the title.

## Architecture

### Wrapper pattern

`getPairedCounts`, `seqSegment`, `getGeneCalls`, and `selectBestMatch`
are symlinks to `execWrapper.sh`. The wrapper exports `SDIR` (repo
directory) and `SNAME` (invoked name) and runs
`bin/rscript.sh $0.R "$@"`. To add a tool: write `foo.R`, symlink
`foo -> execWrapper.sh`, add `docs/usage/foo.doc`.

- `bin/rscript.sh` runs every R script: the wrapped tools, the
  standalone scripts (`fixXChrom.R`, `redactDeliverables.R` from
  `postProcess.sh`; `plotSeg.R` by hand). It puts `Rlib/` first on
  `R_LIBS`, runs `bin/checkRlib.R` against `Rlib.lock.tsv` (exit 1 on
  any mismatch), then `Rscript --vanilla`. Scripts cannot use
  `~/.Rprofile` helpers; `include/tools.R` redefines `cc()` and
  `len()`.
- Wrapped scripts source `include/*.R` via `Sys.getenv("SDIR")`, which
  falls back to the relative path `seqCNA`.
- A new package a script loads goes into `RECORD` in
  `bin/checkRlib.R`, or into `Rlib.lock.tsv` if it computes part of
  the output (see "Pinned versions").

### Shell helpers (bin/)

- `loadTools.sh`: sourced by `PIPE.sh`, `fixChromosomeNames.sh`,
  `postProcess.sh`. Loads `samtools/1.20`, prepends ImageMagick 7,
  checks R 4.5.x, sets `TMPDIR` under `/localscratch/$USER`.
- `rscript.sh`, `checkRlib.R`: the version check (see "Pinned
  versions").
- `getRelease.sh`: prints the release from `VERSION.md`.
  `readme2yaml.R`: `project.yaml` for `delivery.sh`.
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

- One R, 4.5.1, for every script. `seqDNAcopy/` is a git submodule with
  the fork's source (`git submodule update --init`). `00.SETUP.sh`
  builds it into `Rlib/` (gitignored), which holds nothing else, and
  records versions in `Rlib/VERSIONS.txt`. The fork's NAMESPACE imports
  `IRangesList` (0.0.4); `RangesList` is gone from IRanges 2.x.
- pctGCdata 0.3.0 is not on CRAN; the installed copy came from GitHub
  (`veseshan/pctGCdata`, tag v0.3.0, commit d2d4faf).

### Pinned versions

A library update must never change seqCNA output silently.
`Rlib.lock.tsv` (package, version, role) lists only R and the packages
that compute what ends up in the output; the role column says what:
R (`loess` in GC correction, the RNG behind `set.seed`), DNAcopy
(segmentation), pctGCdata (GC content), Rsamtools (`scanBam`, which
reads are counted), data.table (`foverlaps` in getGeneCalls) and
seqDNAcopy. The other packages (dplyr, readr, openxlsx, ...) can make
a run fail but not change its numbers; they are not checked, only
recorded in `Rlib/VERSIONS.txt` (`RECORD` in `bin/checkRlib.R`).
Nothing is installed by the pipeline except seqDNAcopy: the locked
packages come from the shared or personal library.
`bin/rscript.sh` checks the lock before every R script, and `PIPE.sh`
CHECK does once (`rscript.sh --check`); a mismatch stops the run, and
the fix is to reinstall the locked version. `00.RUNSTATUS.txt` records
`R_VERSION`, `RLIB_LOCK_MD5`, `SAMTOOLS_VERSION` and `MAGICK_VERSION`;
samtools and ImageMagick are recorded, not checked.

Versions change only at a planned version scan: install the new
versions, edit `Rlib.lock.tsv` by hand, rerun a reference project,
compare with the previous run, record it in `CHANGELOG.md` and
`docs/methods.md`.

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

### chrX

No chrX copy-number event may be called: that needs a sex-matched
normal and the known sex of the sample, which we almost never have.
`getGeneCalls.R` keeps the autosomes only (1-22 for b37, 1-19 for
mm10), by number, in both the segments and the annotation, so the
rule does not depend on how an annotation names X (Exome `X`,
IMPACT410/505 `23`). X is still counted and segmented: it is in
`rms.logr` (which favors a sex-matched normal in `selectBestMatch`),
the `X.seg.mean.*` QC values, the IGV `.seg` (renamed to X by
`fixXChrom.R`) and the plots. Do not add X to the gene calls.

### Sample ID conventions

`getSampleTag.sh` reads `SM:` from the BAM `@RG` lines and adds an `s_`
prefix if missing. A pair ID is `<tumor>__<normal>`. Internal files
(everything under `out/` and `outAll/`, `bestMatches____out`, `lodir`)
use the pair ID; the deliverables in `<ProjNo>/` other than the PNGs
use the tumor ID. The tumor and normal of a pair are never parsed from
the pair ID in the pipeline: they are columns of `pairs.tsv`, and
`pairing.tsv` maps each best pair to its tumor. A sample can be both a
tumor and a normal, so a regex on sample names is not safe. Pair IDs
are still matched by pattern in `getGeneCalls.R` and
`redactDeliverables.R` (sample columns match `^s_`) and in `plotSeg.R
-r` (`__s_.*$`). Changing the naming breaks them.

### Per-pair outputs

`out/<first 2 chars of tumor ID>/<tumorId>/<pairId>/` contains:
- `<pairId>_Counts.rda` and `.json` (getPairedCounts)
- `<pairId>_seqSeg.rda` (DNAcopy object named `out`), `.seg`, `.png`,
  `.out`, and with `TUMORID=` `.redacted.png` (seqSegment)

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
  the smallest `abs(logR)/sqrt(nProbes)` is set to logR 0. `SEED=42`:
  `set.seed()` before `seqsegment()`, because DNAcopy's p-values are
  by permutation; unseeded, a rerun on the same counts can move
  breakpoints and change gene calls. The seed is in the `.out` file.
- getGeneCalls.R: significant if `FDR < 0.05` and `|seg.mean| > 1`.
  When a gene overlaps several segments, the lowest `seg.mean` is used,
  so deletions take priority over amplifications. For the Exome assay
  the canonical isoforms come from `resources/vcf2maf/v1.6.12/`.
  Autosomes only (see chrX above).
- `docs/methods.md` states these values; update it if they change.
  `docs/output.md` describes every deliverable, its columns, the plot
  and the thresholds; update it if any of them change.

## Environment dependencies

- `~/bin/sbatch`, `~/bin/bsub`, `~/bin/bSync.sh`: the personal Slurm
  wrappers, as in ATAC-seq and ChIP-seq.
  `docs/dev/IRIS_CLUSTER_INFO.md`.
- `module load R/4.5.1_nds` (personal module in `~/.modules`) before
  submitting; batch jobs inherit the environment (`--export=ALL`).
- `samtools/1.20` module; ImageMagick at
  `/usersoftware/core001/common/RHEL_8/imagemagick/7.1.2-33`.
- `/localscratch/$USER` on every node for temp files.
