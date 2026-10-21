# Changelog

All notable changes to this project will be documented in this file.

The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this
project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [6.0.0] - 2026-10-21

The port from JUNO/LSF (CentOS 6/7, R 3.6 and R 4.2) to IRIS/Slurm
(RHEL 8.10, R 4.5.1). Breaking: Slurm only, one R. The specification
and its decisions (D1 to D16) are in
`attic/docs/00.SPEC_IRIS_PORT.md`; how the port was done and validated
is in `docs/dev/SLURM_PORT.md`, with the Proj_12257 run (b37 exome,
6 tumors x 5 normals) it was validated on. The last LSF version is on
branch `juno` and tag `v5-juno` (3e6c18d).

### Breaking
- IRIS/Slurm only. Submit with `sbatch seqCNA/PIPE.sh TARGETTED|sWGS`
  from the project directory; LSF support removed (`bsub`, `bSync`,
  `parseLSF.py`, `-w post_done()`). Logs in `SLURM.*/`. One run per
  project directory; only a run that failed in CHECK can be submitted
  again in the same directory.
- `seqCNA.sh` retired to `attic/`. Its two submissions per pair are the
  COUNT and SEG stages of `PIPE.sh`, which fans each stage out with
  `~/bin/bsub` and waits with `bSync`/`bCheck` between stages.
- One R, 4.5.1 from PATH, for every script. `PIPE.sh` no longer strips
  R from PATH and loads R 3.6.1; `postProcess.sh` no longer loads
  R 4.2.2. The committed R 3.5.1 build of `Rlib/seqDNAcopy` is
  removed.
- Pinned versions of what sets the output. `Rlib.lock.tsv` lists R
  4.5.1, DNAcopy 1.84.0, pctGCdata 0.3.0, Rsamtools 2.26.0,
  data.table 1.18.6.1 and seqDNAcopy 0.0.4, with what each computes.
  Every R script runs through `bin/rscript.sh`, which checks them and
  stops on a mismatch; `PIPE.sh` checks once in CHECK. The other
  packages are not checked, only recorded in `Rlib/VERSIONS.txt`. All
  but seqDNAcopy come from the shared or personal library. The locked
  versions are the ones test20 and the Proj_12257 run used.
- `seqDNAcopy` 0.0.4 (submodule): `NAMESPACE` imports `IRangesList`
  instead of `RangesList`, removed from IRanges 2.x; takes the
  upstream `hg38` addition to `HUMAN.BUILDS`.
- `execWrapper.sh` no longer reads `/etc/system-release`; it runs the
  script through `bin/rscript.sh`.
- Standalone R scripts (`fixXChrom.R`, `redactDeliverables.R`,
  `plotSeg.R`) run through `bin/rscript.sh`, not plain `Rscript`.
- The deliverables name each sample by its tumor ID, not by the pair
  ID `<tumor>__<normal>`, so they do not show which normal was used:
  `<ProjNo>___IGV.seg`, `GeneMatrix.csv`, `SegmentMatrix.csv`,
  `<ProjNo>___GeneTable.xlsx` and `<ProjNo>___seqSeg.pdf`. The file
  names are unchanged. The PNGs in `<ProjNo>/` keep the pair ID and are
  not for delivery. This replaces the manual `redactNormals.R` step,
  which stripped `__.*` by regex and missed the PDF.

### Added
- Run records, from the ChIP-seq port: `00.RUNSTATUS.txt`,
  `SLURM.CTRL/jobs.tsv`, `#SEQCNA_EXIT=<rc>` trailers
  (`bin/runStage.sh`), `bin/checkRun.sh` (exit 0 ok, 1 failed,
  2 running, 3 sacct failed). `COMPLETED` is written only after every
  job is `COMPLETED` in `sacct` and every deliverable is present. A
  failed or cancelled control job marks the run `FAILED` and cancels its
  queued jobs. `00.RUNSTATUS.txt` records `R_VERSION`,
  `RLIB_LOCK_MD5`, `SAMTOOLS_VERSION` and `MAGICK_VERSION`.
- Up-front checks before anything is submitted: BAM directory, every
  BAM readable with an index, `tumors` and `normals` matching the
  relabeled BAM paths, one genome across all BAMs and a supported one,
  project id.
- `out/pipeline_info/pairs.tsv` (PAIR, TUMOR_BAM, NORMAL_BAM, ODIR,
  TUMOR, NORMAL) and `cmd.sh.log` (versions, tools, command line).
- `pairing.tsv` in the run directory (`selectBestMatch`): the best
  normal of each tumor, with TUMOR, NORMAL, PAIR and the `rms.logr` it
  was chosen by. It names the normals and is not delivered.
- `redactDeliverables.R`, run by `postProcess.sh` after `getGeneCalls`:
  replaces each pair ID by its tumor ID in the IGV `.seg`, the two CSVs
  and the xlsx, by the exact map in `pairing.tsv`, not by a regex (a
  sample can be both a tumor and a normal). It stops if a tumor has two
  pairs or an ID is not in the map.
- `seqSegment TUMORID=`: also writes `<pair>_seqSeg.redacted.png`, the
  same plot titled with the tumor ID. `PIPE.sh` sets it, and
  `postProcess.sh` builds the PDF from these plots.
- `seqSegment SEED=` (default 42): `set.seed()` before segmentation,
  written to the `.out` file. DNAcopy computes p-values by
  permutation, and 5.x set no seed: in `../test20` two reruns on the
  same counts changed 3 of 43 best-pair segmentations and 19 gene
  calls. With the seed a rerun is byte-identical (checked on one of
  those pairs; seed 7 instead of 42 gives 52 segments instead of 50).
- Walltime classes SHORT (under two hours, short partitions) and LONG;
  COUNT jobs are sized per pair from the BAM sizes (`runClass`).
- `bin/loadTools.sh`: samtools 1.20 module, ImageMagick 7 from the
  core001 shared tree, R 4.5 check, `TMPDIR` under `/localscratch`.
- `resources/vcf2maf/v1.6.12/`: the two isoform override files
  `getGeneCalls.R` reads for the Exome assay, copied unchanged from
  https://github.com/mskcc/vcf2maf/tree/v1.6.12/data on 2026-10-20
  (they lived under `/opt/common` on JUNO):
  `isoform_overrides_uniprot` (md5 537a33fd15869d1528531e6289f8493e),
  `isoform_overrides_at_mskcc` (md5 cff4b7fc3d1257ad8dd054e0e00101e3).
- **`delivery.sh` delivers on IRIS** (`00.ISSUES.md` #4), as
  ATAC-seq's `deliverResults.sh` does: to
  `/data1/core002/res/bic/results/<pi>/<invest>/Proj_<no>/r_NNN/seqcna`
  (`-d` takes the folder from the path of the analysis directory, the
  run folder is the next `r_NNN`), writing as `svc_core001_bot01`, and
  with `-n` for a dry run. It stops unless `bin/checkRun.sh` exits 0
  and every deliverable is present, and if a delivered table still
  contains `__<normal>` for a normal in `pairing.tsv`. It delivers the
  same files as 5.x (the IGV `.seg`, the PDF, `GeneMatrix.csv`,
  `SegmentMatrix.csv`, the GeneTable xlsx, or `README.txt` when nothing
  is significant) plus `docs/methods.md` and `docs/output.md`; not the
  PNGs or `pairing.tsv`.
  The new `bin/readme2yaml.R` writes `project.yaml` (pipeline `seqcna`,
  the assembly from `GENOME_TAG`) from the first `README.txt` going up
  to `Proj_<no>`, and the run is registered with the bicdelivery
  permissions database. The email (`tmpldeliveryEmail_*.txt`, written
  to `deliveryEmail_<no>_r_NNN_<yymmdd>.txt`) gives the release and run
  folder, links to `docs/output.md` and `docs/methods.md` on GitHub at
  the release tag, and no longer says the BAMs are in the alignment
  folder. The 5.x script copied to `/ifs/res` by hand-given path, with
  no checks. Tested with a delivery of the Proj_12257 run; the normal
  check also on unredacted copies of `GeneMatrix.csv` and the xlsx, and
  the missing-`README.txt` stop on a test directory.
- **`VERSION.md`** holds the client-facing release (`v6.0.0`), set by
  hand when a release is made. `PIPE.sh` records it as `RELEASE=` in
  `00.RUNSTATUS.txt` and `cmd.sh.log` and stops if it is missing
  (`bin/getRelease.sh`); `VERSION=` stays the `git describe` of the
  exact code.
- `docs/dev/SLURM_PORT.md`, `docs/dev/IRIS_CLUSTER_INFO.md`,
  `00.ISSUES.md`, this file.

### Changed
- Segmentation differs from 5.x. DNAcopy 1.75.5 (2023-10-08) fixed
  `segment()` passing the whole weight vector to each chromosome, so
  that every chromosome after chr1 was segmented with chr1's weights.
  `seqsegment()` always passes weights, so every 5.x run, made with an
  older DNAcopy, carries the bug; 6.0.0 uses DNAcopy 1.84.0. On the 38
  tumors x 20 normals of `../test20` the JUNO run made 1,817 gene calls
  and IRIS 1,229 on the same pairs (1,128 shared, none in the opposite
  direction); the IRIS code with DNAcopy 1.60.0 reproduces all 1,817.
  Counts, GC correction and `getGeneCalls` match JUNO. Rerun a 5.x
  project before comparing or merging its calls with 6.x calls. Report:
  `../test20/seqCNA_JUNO_vs_IRIS_concordance.html`.
- `VERSION` in `getPairedCounts.R`, `seqSegment.R` and
  `getGeneCalls.R`, and `Version:` in `DESCRIPTION`, are 6.0.0 (were
  5.dev). `VERSION` is written to the counts `.rda` and the `.out`
  file.
- `PIPE.sh` keeps the LSF first step: run directly with no `tumors`
  or `normals`, it writes `samples` (one line per BAM) and exits 0.
  It now stops before any run record is written, so the run can be
  submitted in the same directory. The early port checked `tumors` and
  `normals` first and so never wrote `samples` in that case.
- `tumors` and `normals` are checked against the relabeled BAM paths,
  as they are applied after FIX. The early port checked the original
  paths, which rejected lines copied from `samples` (`s_<id>__`).
- A sample may be listed in both `tumors` and `normals`; its self-pair
  is skipped. The LSF version ran it and `selectBestMatch` dropped it;
  the early port stopped the run.
- A run that failed in CHECK submitted nothing and no longer blocks the
  directory; a new run replaces its records. The early port refused
  any rerun, so a typo in `tumors` needed a new directory.
- `fixChromosomeNames.sh`: samtools from `bin/loadTools.sh` instead of
  `module load samtools/1.19.2`; calls
  `$SDIR/GenomeData/getGenomeBuildBAM.sh` instead of relying on PATH.
- `GenomeData/getGenomeBuildBAM.sh`: samtools from PATH; exits 1 when
  the header cannot be read instead of printing `unknown` for the md5
  of empty input.
- `postProcess.sh`: `set -eu`; the PNG to PDF merge (`magick`) moved
  here from `PIPE.sh`; `mkdir -p`.
- Resource requests from measurement (2026-10-20, 9 GB exome BAM pairs):
  COUNT and SEG `-n 2 -M 8G` (peak 2.0 and 2.7 GB), FIX `-n 1 -M 4G`,
  POST `-n 1 -M 8G`. Under LSF the pair jobs asked `rusage[mem=16]`
  per slot, advisory; `-M` is a hard cap here.
- `selectBestMatch.R` takes each pair's tumor and normal from
  `pairs.tsv` instead of parsing the pair ID (`__s_`), and keeps exactly
  one normal per tumor (`slice_min`, a tie goes to the first pair ID;
  `top_n` kept every tied pair). `PIPE.sh` requires one best pair per
  tumor, not at least as many as there are tumors.
- `selectBestMatch.R` no longer loads `readxl`; `getGeneCalls.R` no
  longer loads `biomaRt` (neither was used). `fixXChrom.R` and
  `redactDeliverables.R` load the packages they use instead of
  `tidyverse`, which keeps the pinned set small.
- `README.md`, `CLAUDE.md` describe the Slurm pipeline.
- `docs/methods.md` rewritten against the code. It said 10 normals
  were first chosen by insert size; every normal is run. It now
  states the read filters, the GC correction, the sWGS bin size, the
  clustering, the gene annotations and the versions of R, DNAcopy,
  seqDNAcopy, pctGCdata, Rsamtools and data.table.
- New `docs/output.md`: each deliverable, its columns, how to read the
  plot, what the log2 ratios are relative to (no purity or ploidy is
  estimated, so they are not copy numbers) and the calling rule. It is
  delivered with `docs/methods.md`, and the email refers to both.
- `docs/` holds only user-facing docs: `METHODS.md` moved to
  `docs/methods.md` (delivered and linked under that name), the tool
  usage files to `docs/usage/`, the FACETS comparison to
  `docs/facets/`, and the developer docs (`SLURM_PORT.md`,
  `IRIS_CLUSTER_INFO.md`) to `docs/dev/`.
- The FACETS comparison corrected against the code: seqCNA reports log2
  ratios relative to a baseline level of the tumor, with a p-value and
  FDR per segment and no gain or deletion categories; the methods text
  there is replaced by a pointer to `docs/methods.md`; FACETS is cited
  as Shen and Seshan 2016.
- `DELIVERY_ZIP.md` follows the email templates (file names, output.md
  and methods.md, autosomes only).
- Usage text: `seqSegment` `MINBINCOUNT=` says that bins with that many
  normal counts or fewer are removed; `getGeneCalls` `ASSAY=` lists the
  assays. `README.md` copies all three usage texts.

### Fixed
- No gene calls on chrX in any assay. `getGeneCalls.R` keeps the
  autosomes only, by number, in the segments and the annotation. Before,
  X was left out only for the annotations that name it `X` (Exome,
  M-IMPACT_v1) or omit it; IMPACT410 and IMPACT505 name it `23`, as the
  segments do, and so called X genes (AR, ATRX, ... on a JUNO segment
  file of `../test20`). X stays in segmentation, `rms.logr`, the IGV
  `.seg` and the plots. Patch contigs and MT are dropped too.
- `selectBestMatch.R` failed on a project with a single pair: it drops
  the columns that do not vary between pairs, and with one pair that
  included `sampleId`. Under LSF the error was ignored and `out/` was
  not pruned; the port stopped the run at BEST.
- Pruning `out/` after BEST matches `bestMatches____out` on whole lines
  (`grep -vxFf`). The substring match (`fgrep -vf`, also in the LSF
  version) kept the pair `..._s_N10` when the best was `..._s_N1`.
- `PIPE.sh` run directly inside an interactive Slurm job recorded the
  session's `SLURM_JOB_ID` as its control job, so `checkRun.sh`
  reported the run as RUNNING until the session ended. `PIPE.sh` now
  records `RUN_MODE=sbatch` only when it is the batch script itself
  (`$0` is Slurm's `job<id>/slurm_script`); a direct run records
  `CTRL_JOBID=none` and its host and pid, and ignores an inherited
  `SBATCH_SCRIPT_DIR`. `checkRun.sh` checks that the pid is still
  `PIPE.sh`, not only alive. SIGHUP (closed terminal) now marks the run
  FAILED and cancels its jobs, as SIGTERM does.

### Removed
- `testR`, `testR.R`, `docs/testR.doc` (empty), `misc/INSTALL.RLIB`,
  `Rlib/INSTALL.md`, `GenomeData/genomeInfo_*.sh` (all `/ifs` paths,
  never sourced), to `attic/`.
- `redactNormals.R`, to `attic/`: the pipeline writes redacted
  deliverables.
- `DELIVERY.md`: the email is now built by `delivery.sh` from
  `tmpldeliveryEmail_*.txt`. `DELIVERY_ZIP.md` stays, for a delivery by
  email.

### Not tested
- The mouse path (`mm10` relabel, `Mouse-Twist` assay) and sWGS
  (`BINSIZE=auto`) were not run on IRIS; no such project was at hand
  (`00.ISSUES.md` #2).
- The `delivery.sh` email for a run with no calls (`README.txt`
  instead of the tables) has not been seen on a run.
- `RELEASE=` in `00.RUNSTATUS.txt`: no run has been made with it; the
  Proj_12257 run predates it.
