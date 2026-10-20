# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

To be released as 6.0.0: the port from JUNO/LSF (CentOS 6/7, R 3.6 and
R 4.2) to IRIS/Slurm (RHEL 8.10, R 4.5.1). Breaking: Slurm only, one R.
The specification and its decisions (D1 to D16) are in
`00.SPEC_IRIS_PORT.md`; how the port was done and validated is in
`docs/SLURM_PORT.md`. The last LSF version is tagged `v5-juno`.

### Breaking
- IRIS/Slurm only. Submit with `sbatch seqCNA/PIPE.sh TARGETTED|sWGS`
  from the project directory; LSF support removed (`bsub`, `bSync`,
  `parseLSF.py`, `-w post_done()`). Logs in `SLURM.*/`. One run per
  project directory.
- `seqCNA.sh` retired to `ATTIC/`. Its two submissions per pair are the
  COUNT and SEG stages of `PIPE.sh`, which fans each stage out with
  `~/bin/bsub` and waits with `bSync`/`bCheck` between stages.
- One R, 4.5.x from PATH, for every script. `PIPE.sh` no longer strips
  R from PATH and loads R 3.6.1; `postProcess.sh` no longer loads
  R 4.2.2. The R packages the pipeline needs beyond the shared install
  (DNAcopy, pctGCdata, RJSONIO) come from the user's library;
  `seqDNAcopy` is built once per checkout into `Rlib/` (gitignored) by
  `00.SETUP.sh`. The committed R 3.5.1 build of `Rlib/seqDNAcopy` is
  removed.
- `seqDNAcopy` 0.0.4 (submodule): `NAMESPACE` imports `IRangesList`
  instead of `RangesList`, removed from IRanges 2.x; takes the
  upstream `hg38` addition to `HUMAN.BUILDS`.
- `execWrapper.sh` no longer reads `/etc/system-release`; it requires
  R 4.5.x and `Rlib/seqDNAcopy`.
- `redactNormals.R` calls `openxlsx::write.xlsx` directly instead of
  `write_xlsx()` from `~/.Rprofile`; run standalone R scripts with
  `Rscript --no-init-file`.

### Added
- Run records, from the ChIP-seq port: `00.RUNSTATUS.txt`,
  `SLURM.CTRL/jobs.tsv`, `#SEQCNA_EXIT=<rc>` trailers
  (`bin/runStage.sh`), `bin/checkRun.sh` (exit 0 ok, 1 failed,
  2 running, 3 sacct failed). `COMPLETED` is written only after every
  job is `COMPLETED` in `sacct` and every deliverable is present. A
  failed or cancelled control job marks the run `FAILED` and cancels its
  queued jobs.
- Up-front checks before anything is submitted: BAM directory, every
  BAM readable with an index, `tumors` and `normals` present and
  matching, one genome across all BAMs and a supported one, project id.
- `out/pipeline_info/pairs.tsv` (PAIR, TUMOR_BAM, NORMAL_BAM, ODIR) and
  `cmd.sh.log` (versions, tools, command line).
- Walltime classes SHORT (under two hours, short partitions) and LONG;
  COUNT jobs are sized per pair from the BAM sizes (`runClass`).
- `bin/loadTools.sh`: samtools 1.20 module, ImageMagick 7 from the
  core001 shared tree, R 4.5 check, `TMPDIR` under `/localscratch`.
- `resources/vcf2maf/v1.6.12/`: the two isoform override files
  `getGeneCalls.R` reads for the Exome assay, vendored from vcf2maf
  v1.6.12 (they lived under `/opt/common` on JUNO).
- `docs/SLURM_PORT.md`, `docs/IRIS_CLUSTER_INFO.md`, `00.ISSUES.md`,
  this file.

### Changed
- `fixChromosomeNames.sh`: samtools from `bin/loadTools.sh` instead of
  `module load samtools/1.19.2`; calls `$SDIR/GenomeData/getGenomeBuildBAM.sh`
  instead of relying on PATH.
- `GenomeData/getGenomeBuildBAM.sh`: samtools from PATH; exits 1 when
  the header cannot be read instead of printing `unknown` for the md5
  of empty input.
- `postProcess.sh`: `set -eu`; the PNG to PDF merge (`magick`) moved
  here from `PIPE.sh`; `mkdir -p`.
- Resource requests from measurement (2026-10-20, 9 GB exome BAM pairs):
  COUNT and SEG `-n 2 -M 8G` (peak 2.0 and 2.7 GB), FIX `-n 1 -M 4G`,
  POST `-n 1 -M 8G`. Under LSF the pair jobs asked `rusage[mem=16]`
  per slot, advisory; `-M` is a hard cap here.
- `selectBestMatch.R` no longer loads `readxl`; `getGeneCalls.R` no
  longer loads `biomaRt` (neither was used).
- `README.md`, `CLAUDE.md` describe the Slurm pipeline.

### Removed
- `testR`, `testR.R`, `misc/INSTALL.RLIB`, `Rlib/INSTALL.md`,
  `GenomeData/genomeInfo_*.sh` (all `/ifs` paths, never sourced), to
  `ATTIC/`.

### Not tested
- The mouse path (`mm10` relabel, `Mouse-Twist` assay) and sWGS
  (`BINSIZE=auto`) were not run on IRIS; no such project was at hand
  (`00.ISSUES.md`).
- Output was not compared with a JUNO run of the same BAMs
  (`00.ISSUES.md`).
