# Porting seqCNA from JUNO/LSF to IRIS/Slurm

Done 2026-10-20 on `feat/iris`, for 6.0.0. The survey, design and
decisions (D1 to D16) are in `attic/docs/00.SPEC_IRIS_PORT.md`. The
model is the ChIP-seq port (`docs/SLURM_PORT.md` there, v1.0.0),
itself from ATAC-seq: a synchronous control job that fans each stage
out with `~/bin/bsub` and blocks with `bSync`/`bCheck`. The PEMapper pattern
(the whole DAG submitted at once with `afterok` ids) was considered and
not used (spec 3.3). JUNO is gone; there is no LSF path left
(`attic/seqCNA.sh` is reference only).

## Strategy: keep the pipeline shape, change three things

1. How R is found: one R 4.5.x from PATH for every script, `seqDNAcopy`
   built per checkout into `Rlib/` (`00.SETUP.sh`), the common packages
   from the user's library. No `module load R/...`, no `/etc/system-release`.
2. How jobs are submitted and waited on: the ChIP-seq control-job
   primitives, with `CHIP_` renamed `SEQCNA_`.
3. Where tools and resources live: `bin/loadTools.sh`, `resources/vcf2maf/`.

| LSF | Slurm | Where |
| --- | --- | --- |
| `nohup ./seqCNA/PIPE.sh TARGETTED` (or `bsub`) | `sbatch ./seqCNA/PIPE.sh TARGETTED` | `~/bin/sbatch`, which exports `SBATCH_SCRIPT_DIR` |
| `bsub` | `~/bin/bsub` | SchedMD openlava shim with local patches |
| `bSync NAME` | `bSync.sh "^NAME$"` | `~/bin/bSync.sh`, wrapped in `bin/slurmTools.sh` |
| `parseLSF.py LSF*/* \| fgrep -v Succ` | `bCheck NAME NJOBS` (`sacct` State), `bCheckAll REGEX` | `bin/slurmTools.sh` |
| `bsub -w "post_done(COUNT)" seqSegment` | COUNT stage, `waitStage`, SEG stage | `PIPE.sh` |

## Stages

| Stage | Jobs | Command | Log dir | Class | Resources |
| --- | --- | --- | --- | --- | --- |
| CHECK | control job | inputs, genome, tumors/normals, project id | | | |
| FIX | one per BAM | `fixChromosomeNames.sh BAM` | `SLURM.01.FIX` | SHORT | `-n 1 -M 4G` |
| COUNT | one per tumor x normal | `getPairedCounts GENOME= NORMAL= TUMOR= ODIR= SAMPLEID=` | `SLURM.02.COUNT` | `runClass`, 1 min/GB | `-n 2 -M 8G` |
| SEG | one per pair | `seqSegment BINSIZE= ODIR= COUNTS=` | `SLURM.03.SEG` | SHORT | `-n 2 -M 8G` |
| BEST | control job | `selectBestMatch out`, `outAll/` hard links, prune `out/` | | | |
| POST | one | `postProcess.sh ASSAY PROJNO` (incl. the PDF) | `SLURM.04.POST` | SHORT | `-n 1 -M 8G` |
| DELIVERABLES, STAGING | control job | file checks, `bCheckAll`, `COMPLETED` | | | |

Job names are `qseqCNA_<STAGE>_<pid>`; every `sacct` query is bounded
by `SEQCNA_START`. The control job: `cmobic_cpu`, `--qos=priority`,
3 days, `-c 2 --mem=16G`, log `SLURM.CTRL/<jobid>.out`.

## Flag translation

| LSF (v5, `seqCNA.sh` / `PIPE.sh`) | Slurm (through the shim) | Note |
| --- | --- | --- |
| FIX `-W 24:00` | `-n 1 -M 4G`, SHORT | LSF `24:00` is 24 hours; to the shim it would be 24 minutes. The job links or reheaders one BAM: seconds |
| COUNT `-W 24:00 -n 2 -R "rusage[mem=16]"` | `-n 2 -M 8G`, `runClass` | measured: 176 s and 2.0 GB peak RSS for a 9.0 GB exome pair |
| SEG `-W 24:00 -n 2 -R "rusage[mem=16]" -w "post_done(...)"` | `-n 2 -M 8G`, SHORT | measured: 66 s, 2.7 GB. `-w` removed: Slurm needs numeric ids and an `afterok` on a failed parent pends forever |
| `convert` in `PIPE.sh` | `magick` in `postProcess.sh` | ImageMagick 7 from `/usersoftware/core001/common/RHEL_8/imagemagick/7.1.2-33` |

`-M` is the job's total memory and a hard cgroup cap; LSF `rusage` was
per slot and advisory.

## Walltime classes

As ChIP-seq: SHORT (`-W 1:59:00`, `cmobic_short,cpushort`, no qos) or
LONG (`-W 12:00:00`, `cmobic_cpu`, qos priority on that one `bsub`
call). COUNT picks per pair with `runClass MIN_PER_GB TUMOR NORMAL`:
SHORT if the estimate is at most 60 min (`SEQCNA_SHORT_MAX_MIN`).

| Rate | min/GB | Input | Measured |
| --- | --- | --- | --- |
| `RATE_COUNT` | 1 | tumor + normal BAM | 0.33 min/GB in the foreground on a login node, 2026-10-20 (176 s, 9.0 GB); 3x margin |

A 9 to 10 GB exome pair estimates at 9 to 10 min, SHORT. An sWGS pair
would have to exceed 60 GB before it went LONG.

## Run records and checking a run

As ChIP-seq, with `SEQCNA` in place of `CHIP`: `00.RUNSTATUS.txt`
(adds `ASSAY_TYPE`, `BINSIZE`, `BAMDIR`, `GENOME_TAG`, `GENOME`, `ASSAY`,
`PROJNO`, `NPAIRS`, `SAMPLES`), `SLURM.CTRL/jobs.tsv`,
`#SEQCNA_EXIT=<rc>` trailers, `bin/checkRun.sh` (exit 0 ok, 1 failed,
2 running, 3 sacct failed). `MESSAGE` holds the reason for a failure:
the failed up-front check, or the first FATAL line of each failed stage
job. `STATUS=COMPLETED` is written last, after the `sacct` sweep and the
deliverables check. One run per project directory: `out`, `outAll`,
`bamRelabel`, `LSF`, `LSF.FIX` or the run records block a new run. The
exception is a run that failed in CHECK: it submitted nothing, and a
new run replaces its records.

```bash
./seqCNA/bin/checkRun.sh
grep -L "#SEQCNA_EXIT=0" SLURM.0*/*.out
```

## Environment changes

| What | Was | Now |
| --- | --- | --- |
| R for counting, segmentation, best match | `module load R/R-3.6.1`, `R_LIBS=/home/socci/lib/R/CentOS7/3.6.1`, R stripped from PATH first | R 4.5.1 from PATH (`module load R/4.5.1_nds`), `Rscript --vanilla` with `R_LIBS=$SDIR/Rlib:...` |
| R for `fixXChrom.R`, `getGeneCalls` | `module load R/R-4.2.2` in `postProcess.sh` | the same R 4.5.1; `Rscript --no-init-file` |
| `seqDNAcopy` | committed `Rlib/seqDNAcopy` built with R 3.5.1 | built by `00.SETUP.sh` from the submodule (0.0.4, NAMESPACE fix) into gitignored `Rlib/` |
| DNAcopy, pctGCdata, RJSONIO | `/home/socci/lib/R/...` | user's R library (`~/R/x86_64-pc-linux-gnu-library/4.5`); `00.SETUP.sh` checks them |
| samtools | `module load samtools/1.19.2`; `/opt/common/.../samtools-1.9` hard-coded | `module load samtools/1.20` in `bin/loadTools.sh` |
| ImageMagick | `convert` on PATH | `magick` from the core001 shared tree (`bin/loadTools.sh`) |
| vcf2maf isoform overrides | `/opt/common/CentOS_6-dev/vcf2maf/v1.6.12/data/` | `resources/vcf2maf/v1.6.12/` |
| `write_xlsx()` | `~/.Rprofile` | `openxlsx::write.xlsx` |
| TMPDIR | node `/tmp` | per job under `${SEQCNA_LOCAL_TMP:-/localscratch/$USER}` (`bin/runStage.sh`) |
| OS check | `/etc/system-release` major 6 or 7 in `execWrapper.sh` | none |

## Validation

Done 2026-10-20 on `feat/iris` (uncommitted, so the runs report
`proj/B-100-635-54-g06971b0-dirty`), against
`attic/docs/00.SPEC_IRIS_PORT.md` section 4. Inputs: the Proj_12257
BAMs (b37 exome, 5 tumors x 2 normals, 4.0 to 5.1 GB each, full size;
no downsampling was needed at these run times). Test directories are beside the checkout in `../test01/`,
`../test02/`, `../test03/`.

| Item | Test | Result |
| --- | --- | --- |
| 4.1 setup | `00.SETUP.sh` | seqDNAcopy 0.0.4 built and loaded under R 4.5.1; `Rlib/VERSIONS.txt` |
| 4.2 hand run | `getPairedCounts` then `seqSegment` on one pair, foreground | 176 s / 2.0 GB and 66 s / 2.7 GB; `.seg`, `.out`, `.png` as before (60,292 bins, 27 segments) |
| 4.3, 4.4 full run | test01, control job 546394 | `COMPLETED`, 28 of 28 jobs `COMPLETED` (7 FIX, 10 COUNT, 10 SEG, 1 POST), every log `#SEQCNA_EXIT=0`, `checkRun.sh` 0; 6 min 10 s wall, all stage jobs on `cmobic_short`. Deliverables: 5 PNGs, `___IGV.seg` (180 rows, X renamed), `GeneMatrix.csv` (8 genes), `SegmentMatrix.csv`, `___GeneTable.xlsx`, `___seqSeg.pdf` (5 pages). `out/` pruned to the 5 best pairs, `outAll/` holds all 10 |
| failure: missing index | test02, one BAM without `.bai`, `PIPE.sh` run directly | `FAILED` at `STAGE=CHECK`, rc 1, no job submitted, `MESSAGE=no index (.bai) for [...]` |
| failure: rerun | test01, `PIPE.sh` again | refused before the EXIT trap: "this directory already holds a run: 00.RUNSTATUS.txt SLURM.CTRL/jobs.tsv out outAll bamRelabel"; the first run's records untouched |
| failure: cancelled stage job | test03, control job 588179, `scancel` of COUNT job 588659 while running | `FAILED` at `STAGE=COUNT`, rc 1, after the other 9 COUNT jobs finished; `MESSAGE=COUNT failed: <pair> job 588659 CANCELLED: ... DUE to SIGNAL Terminated`; no SEG job submitted; `checkRun.sh` 1 with the job listed under cancelled |
| redactNormals | copy of the test01 deliverables, `Rscript --no-init-file` | rc 0; `__<normal>` removed from the four tables, xlsx rewritten by `openxlsx::write.xlsx` |
| 4.5 comparison with JUNO | | not done: Proj_12257 has no JUNO run (`00.ISSUES.md` #1) |
| 4.6 mouse, sWGS | | not done (`00.ISSUES.md` #2) |

sacct for test01 (`MaxRSS` is the cgroup watermark and includes page
cache, see `docs/IRIS_CLUSTER_INFO.md`):

| Stage | Jobs | Elapsed | MaxRSS | ReqMem |
| --- | --- | --- | --- | --- |
| FIX | 546419 to 546430 | 2 s | 47 to 58 MB | 4G |
| COUNT | 546510 to 546532 | 2:52 to 3:24 | 1.9 to 8.4 GB (cache-inflated; 2.0 GB real in the foreground) | 8G |
| SEG | 548253 to 548272 | 1:00 to 1:22 | 2.3 to 2.9 GB | 8G |
| POST | 561199 | 9 s | 459 MB | 8G |
| control | 546394 | 6:10 | 190 MB | 16G |
