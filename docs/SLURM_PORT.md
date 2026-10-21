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

1. How R is found: one R 4.5.1 from PATH for every script,
   `seqDNAcopy` from `Rlib/` (built per checkout by `00.SETUP.sh`),
   the other packages from the shared or personal library. R and the
   packages that set the output are checked against `Rlib.lock.tsv`
   before every R script (`bin/rscript.sh`). No `module load R/...`,
   no `/etc/system-release`.
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
| `RATE_COUNT` | 1 | tumor + normal BAM | 0.33 min/GB, 2026-10-20: in the foreground on a login node (176 s, 9.0 GB) and on every Proj_12257 COUNT job (4:39 for 14.3 GB to 9:15 for 28.3 GB); 3x margin |

A 9 to 10 GB exome pair estimates at 9 to 10 min, SHORT. An sWGS pair
would have to exceed 60 GB before it went LONG.

## Run records and checking a run

As ChIP-seq, with `SEQCNA` in place of `CHIP`: `00.RUNSTATUS.txt`
(adds `ASSAY_TYPE`, `BINSIZE`, `BAMDIR`, `GENOME_TAG`, `GENOME`, `ASSAY`,
`PROJNO`, `NPAIRS`, `SAMPLES`), `SLURM.CTRL/jobs.tsv`,
`#SEQCNA_EXIT=<rc>` trailers, `bin/checkRun.sh` (exit 0 ok, 1 failed,
2 running, 3 sacct failed). `RUN_MODE` is `sbatch` when `PIPE.sh` is
the batch script of `$SLURM_JOB_ID` (`$0` is
`/var/spool/slurmd/job<id>/slurm_script`, checked with job 1264967),
and `direct` otherwise, including a run started by hand inside an
interactive job: there `SLURM_JOB_ID` names the session, and the early
port recorded it as `CTRL_JOBID`, so `checkRun.sh` reported the run as
RUNNING for as long as the session lasted. A direct run records
`CTRL_JOBID=none`; `checkRun.sh` checks that `PID` on `HOST` is still
`PIPE.sh`. `MESSAGE` holds the reason for a failure:
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
| R for counting, segmentation, best match | `module load R/R-3.6.1`, `R_LIBS=/home/socci/lib/R/CentOS7/3.6.1`, R stripped from PATH first | R 4.5.1 from PATH (`module load R/4.5.1_nds`), `bin/rscript.sh`: `Rscript --vanilla` with `Rlib/` first on `R_LIBS`, checked against `Rlib.lock.tsv` |
| R for `fixXChrom.R`, `getGeneCalls` | `module load R/R-4.2.2` in `postProcess.sh` | the same R 4.5.1 through `bin/rscript.sh` |
| `seqDNAcopy` | committed `Rlib/seqDNAcopy` built with R 3.5.1 | built by `00.SETUP.sh` from the submodule (0.0.4, NAMESPACE fix) into gitignored `Rlib/` |
| DNAcopy, pctGCdata, RJSONIO and the rest | `/home/socci/lib/R/...` | the shared or personal library; R, DNAcopy, pctGCdata, Rsamtools and data.table checked against `Rlib.lock.tsv` |
| samtools | `module load samtools/1.19.2`; `/opt/common/.../samtools-1.9` hard-coded | `module load samtools/1.20` in `bin/loadTools.sh` |
| ImageMagick | `convert` on PATH | `magick` from the core001 shared tree (`bin/loadTools.sh`) |
| vcf2maf isoform overrides | `/opt/common/CentOS_6-dev/vcf2maf/v1.6.12/data/` | `resources/vcf2maf/v1.6.12/` |
| `write_xlsx()` | `~/.Rprofile` | `openxlsx::write.xlsx` |
| TMPDIR | node `/tmp` | per job under `${SEQCNA_LOCAL_TMP:-/localscratch/$USER}` (`bin/runStage.sh`) |
| OS check | `/etc/system-release` major 6 or 7 in `execWrapper.sh` | none |

## Validation

Done 2026-10-20 against `attic/docs/00.SPEC_IRIS_PORT.md` section 4,
in two parts.

The Proj_12257 run is in `../../CNA/` (beside `SeqCNAPort/`): the
project's 11 BAMs (b37 exome, 6 tumors and 5 normals, 5.4 to 15.6 GB
each), run as a new project is run: `PIPE.sh` first without `tumors`
and `normals`, which writes `samples`, then the two files built from
its lines. The clone was `devs` f231abf, the 6.0.0 code before the
`VERSION` bump, with every package but seqDNAcopy from the personal
library.

The test runs exercise the pipeline, not the calls. Their input is
`SeqCNAPort/pipeline/bam`: seven normal BAMs (`s_C_*_N001_d01`, 4.0 to
5.1 GB, full size) linked from the argos output under
`/data1/core002/res/ci/argos/18473/`, listed as tumors as well as
normals. test01 to test03 ran on `feat/iris` uncommitted (so they
report `proj/B-100-635-54-g06971b0-dirty`) with DNAcopy and pctGCdata
from temporary copies in `Rlib/`; test10 to test14 ran from fresh
`devs` clones with the personal library. The test directories were
beside the checkout (`../test01/` and so on) and were deleted on
2026-10-20; the job ids can still be checked with `sacct`.

| Item | Test | Result |
| --- | --- | --- |
| 4.3, 4.4 Proj_12257 | `../../CNA`, run directly, pid 3316002 on `isca006` | `COMPLETED`, 72 of 72 jobs (11 FIX, 30 COUNT, 30 SEG, 1 POST), every log `#SEQCNA_EXIT=0`, `checkRun.sh` 0; 11 min 35 s wall, all stage jobs on `cmobic_short`. Best normal `s_C_XCC0PR_N001_d` for all 6 tumors, the lowest `rms.logr` in each case (checked against the 30 `.out` files). Deliverables: 6 PNGs, `___IGV.seg` (306 segments), `GeneMatrix.csv` (21 genes, all deletions, in 2 tumors), `SegmentMatrix.csv`, `___GeneTable.xlsx`, `___seqSeg.pdf`. `out/` pruned to the 6 best pairs, `outAll/` holds all 30 |
| 4.1 setup | `00.SETUP.sh` | seqDNAcopy 0.0.4 built and loaded under R 4.5.1; `Rlib/VERSIONS.txt` |
| 4.2 hand run | `getPairedCounts` then `seqSegment` on one pair of test BAMs, foreground | 176 s / 2.0 GB and 66 s / 2.7 GB; `.seg`, `.out`, `.png` as before (60,292 bins, 27 segments) |
| full run, test BAMs | test01 (5 as tumors x 2 normals), control job 546394 | `COMPLETED`, 28 of 28 jobs `COMPLETED` (7 FIX, 10 COUNT, 10 SEG, 1 POST), every log `#SEQCNA_EXIT=0`, `checkRun.sh` 0; 6 min 10 s wall, all stage jobs on `cmobic_short`. Deliverables: 5 PNGs, `___IGV.seg` (180 rows, X renamed), `GeneMatrix.csv` (8 genes), `SegmentMatrix.csv`, `___GeneTable.xlsx`, `___seqSeg.pdf` (5 pages). `out/` pruned to the 5 best pairs, `outAll/` holds all 10 |
| failure: missing index | test02, one BAM without `.bai`, `PIPE.sh` run directly | `FAILED` at `STAGE=CHECK`, rc 1, no job submitted, `MESSAGE=no index (.bai) for [...]` |
| failure: rerun | test01, `PIPE.sh` again | refused before the EXIT trap: "this directory already holds a run: 00.RUNSTATUS.txt SLURM.CTRL/jobs.tsv out outAll bamRelabel"; the first run's records untouched |
| failure: cancelled stage job | test03, control job 588179, `scancel` of COUNT job 588659 while running | `FAILED` at `STAGE=COUNT`, rc 1, after the other 9 COUNT jobs finished; `MESSAGE=COUNT failed: <pair> job 588659 CANCELLED: ... DUE to SIGNAL Terminated`; no SEG job submitted; `checkRun.sh` 1 with the job listed under cancelled |
| redactNormals | copy of the test01 deliverables, `Rscript --no-init-file` | rc 0; `__<normal>` removed from the four tables, xlsx rewritten by `openxlsx::write.xlsx` |
| one tumor, one normal | test11 (`devs` 95deb1c), control job 1236320 | `COMPLETED`, 10 of 10 jobs (7 FIX, 1 COUNT, 1 SEG, 1 POST), `checkRun.sh` 0, 5 min 41 s. `bestMatches____out` holds the one pair, `out/` not pruned; `README.txt` ("No significant copy number changes") instead of the gene tables. The `.seg` is byte-identical to the same pair in test01 |
| self-pairs skipped | test10 (`devs` 95deb1c), seven samples in both `tumors` and `normals`, run directly inside interactive job 168961 | `COMPLETED`, 92 jobs (7 FIX, 42 COUNT, 42 SEG, 1 POST); no self-pair in `pairs.tsv`, `out/` (7 best) or `outAll/` (42). Run before the `RUN_MODE` fix, so it recorded `CTRL_JOBID=168961` |
| direct run | test12, `PIPE.sh` run directly, pid 3105705 on `isca006` | `COMPLETED`, 10 of 10 jobs, `RUN_MODE=direct`, `CTRL_JOBID=none`, `checkRun.sh` 0; `___IGV.seg` identical to test11 |
| failure: SIGHUP | test13, direct run, SIGHUP sent during COUNT | `FAILED` at `STAGE=COUNT`, rc 129, COUNT job 1265961 cancelled; `checkRun.sh` 2 while running, 1 after |
| sbatch run | test14, control job 1265859 | `COMPLETED`, 10 of 10 jobs, `RUN_MODE=sbatch`, `checkRun.sh` 0; `___IGV.seg` identical to test11 |
| 4.5 comparison with JUNO | `../test20` (`devs` f231abf), a JUNO project rerun: 38 tumors x 20 normals, 760 pairs, run directly, pid 3122405 on `isca006` | `COMPLETED`, 1,579 of 1,579 jobs (58 FIX, 760 COUNT, 760 SEG, 1 POST), 19 min 50 s wall. Counts, bin grid and `getGeneCalls` match JUNO exactly; segmentation does not (gene calls on the same pairs: 1,817 JUNO, 1,229 IRIS, 1,128 shared, none opposite), because DNAcopy 1.75.5 fixed the weights in weighted segmentation; DNAcopy 1.60.0 reproduces every JUNO call. Closed as explained (`00.ISSUES.md` #1); report `../test20/seqCNA_JUNO_vs_IRIS_concordance.html` |
| 4.6 mouse, sWGS | | not done (`00.ISSUES.md` #2) |
| pinned library, seed | `00.SETUP.sh` (60 packages from source, about 23 min on `isca006`), then test31 and test32 (`feat/pin-versions` uncommitted on 2ce47ac, test30's inputs: 5 tumors x 3 normals), control jobs 1325450 and 1325451 | `COMPLETED`, 36 of 36 jobs each, `checkRun.sh` 0; `00.RUNSTATUS.txt` records `R_VERSION=4.5.1`, `RLIB_LOCK_MD5`, `SAMTOOLS_VERSION=1.20`, `MAGICK_VERSION=7.1.2-33`. Every file under `outAll/` (counts, `.rda`, `.seg`, `.out`, PNGs), `pairing.tsv`, `GeneMatrix.csv`, `SegmentMatrix.csv`, xlsx and PDF identical between the two runs; `___IGV.seg` has the same rows in a different order (`find` order in `postProcess.sh`). `bin/checkRlib.R` stops on a wrong R or package version and on a package loaded from another library |
| reduced lock | `Rlib.lock.tsv` cut to R, DNAcopy, pctGCdata, Rsamtools, data.table and seqDNAcopy, the other packages from the shared or personal library; `00.SETUP.sh` builds only seqDNAcopy (seconds). test33 (test30's inputs), control job 1528228 | `COMPLETED`, 36 of 36 jobs, about 7 min wall. Identical to test31 except the git tag (`.out`, `.json`, `GITTag` in `_Counts.rda`) and the timestamps inside the xlsx and PDF (xlsx cells identical). `bin/checkRlib.R` stops on a wrong locked version and on anything but seqDNAcopy in `Rlib/` |

sacct for Proj_12257 (`MaxRSS` is the cgroup watermark and includes
page cache, see `docs/IRIS_CLUSTER_INFO.md`; a direct run has no
control job):

| Stage | Jobs | Elapsed | MaxRSS | ReqMem |
| --- | --- | --- | --- | --- |
| FIX | 1278043 to 1278056 | 1 to 2 s (b37: links only) | 46 to 61 MB | 4G |
| COUNT | 1278119 to 1278157 | 4:39 to 9:15 | 1.8 to 8.4 GB (22 of 30 at the cap; cache-inflated) | 8G |
| SEG | 1280450 to 1280491 | 0:44 to 1:08 | 1.9 to 2.9 GB | 8G |
| POST | 1280960 | 5 s | 359 MB | 8G |

sacct for test20 (BAMs 0.5 to 2.3 GB; the genome tag `b37+mm10` is
not `b37`, so FIX reheaders instead of linking):

| Stage | Jobs | Elapsed | MaxRSS | ReqMem |
| --- | --- | --- | --- | --- |
| FIX | 1269562 to 1269675 | 0:15 to 1:06 | 1.0 to 4.0 GB (4 of 58 at the cap; cache-inflated) | 4G |
| COUNT | 1270334 to 1272829 | 0:41 to 1:49 | 1.3 to 5.3 GB | 8G |
| SEG | 1273218 to 1275214 | 0:18 to 1:12 | 0.8 to 2.7 GB | 8G |
| POST | 1275364 | 10 s | 1.3 GB | 8G |

sacct for test01:

| Stage | Jobs | Elapsed | MaxRSS | ReqMem |
| --- | --- | --- | --- | --- |
| FIX | 546419 to 546430 | 2 s | 47 to 58 MB | 4G |
| COUNT | 546510 to 546532 | 2:52 to 3:24 | 1.9 to 8.4 GB (cache-inflated; 2.0 GB real in the foreground) | 8G |
| SEG | 548253 to 548272 | 1:00 to 1:22 | 2.3 to 2.9 GB | 8G |
| POST | 561199 | 9 s | 459 MB | 8G |
| control | 546394 | 6:10 | 190 MB | 16G |
