# IRIS cluster facts

Measured, not assumed. Probed on IRIS (MSKCC) on 2026-10-19 and
2026-10-20 from `isca006` during the seqCNA Slurm port, on top of the
PEMapper survey of 2026-09 (`PEMapper/docs/IRIS_CLUSTER_INFO.md`, which
has the scheduler settings and the `sacct` state tests in full). Job
ids are kept so a claim can be re-checked with `sacct`.

Design decisions that follow from these facts are in
`attic/docs/00.SPEC_IRIS_PORT.md`; open work in `00.ISSUES.md`. This
file is facts only.

## Scheduler

Slurm, RHEL 8.10. Account `core001`.

| Setting | Value | Consequence |
|---|---|---|
| Default partition | `cpu` | Denies `core001`; `-p` must be explicit. `~/bin/bsub` picks `cmobic_short,cpushort` under two hours, `cmobic_cpu` otherwise. |
| `--mem` | hard cgroup cap | LSF `rusage[mem=]` was advisory. |
| `EnforcePartLimits` | `ALL` | A SHORT job with a qos is rejected on `cpushort`. |
| `MaxArraySize` | 10001 | Not used by this pipeline. |

Slurm writes nothing of its own into the `-o` file. "Did it work" is
`sacct` State (`COMPLETED`), never ExitCode (OOM is `OUT_OF_MEMORY`
with `ExitCode=0:125`).

| Partition | `MaxTime` | Notes |
|---|---|---|
| `cpushort` | `02:00:00` | Pins `cpushort_qos`; no `--qos`. |
| `cmobic_short` | `03:00:00` | |
| `cmobic_cpu` | `7-00:00:00` | Control job, with `--qos=priority`. |

## Personal wrappers in ~/bin

The pipeline depends on them, as ATAC-seq and ChIP-seq do.

| Wrapper | What it does |
|---|---|
| `~/bin/sbatch` | Writes a temp copy of the batch script with `export SBATCH_SCRIPT_DIR=<dir of the script>` inserted after the `#SBATCH` header, then calls `/usr/bin/sbatch`. `PIPE.sh` takes `SDIR` from it. Resolves symlinks (`readlink -f`), so `seqCNA -> ../seqCNA` works. |
| `~/bin/bsub` | SchedMD LSF-to-Slurm shim, patched: `-W` to `-t` (sbatch time format, so `-W 24:00` is 24 minutes), `-M` to `--mem`, `-n` to `-c`, `-o DIR/` to `DIR/%j.out`, partition from the run time. `-w` maps to `afterok` and needs numeric ids; not used. |
| `~/bin/bSync.sh REGEX` | Polls `squeue` until no job of the user matches the name regex. |

## Filesystems

- `/data1/core001/work/bic/socci` (`~/Work`): project and checkout
  storage, visible from every node.
- `/localscratch`: node-local disk (2.8 TB on `isca006`); per-job
  `TMPDIR` goes here (`bin/runStage.sh`). Not shared: a run directory
  or an R library here is invisible to the compute nodes.
- `/tmp`: a 137 GB volume shared by every job on a node. Nothing of
  this pipeline writes there.

## Toolchain

| Tool | Version | Where |
|---|---|---|
| R | 4.5.1 | `/usersoftware/core001/common/RHEL_8/gcc-12.2.0/R/4.5.1`, personal module `R/4.5.1_nds` in `~/.modules`. Sets `R_LIBS_USER=~/R/x86_64-pc-linux-gnu-library/4.5`, which is `/data1/core001/work/bic/socci/opt/lib/R/...`. |
| R packages in the shared or personal library | tidyverse 2.0.0, Rsamtools 2.26.0, IRanges 2.44.0, data.table 1.18.6.1, openxlsx 4.2.9, jsonlite 2.0.0, Rcpp 1.1.2, RJSONIO 2.0.5 | Missing from both as of 2026-10-20: DNAcopy, pctGCdata (user to install). |
| Other R modules | `r/3.6.3`, `R/4.3.0` (empty libraries), `r/4.4.1` (partial) | Not used. |
| samtools | 1.20 | module `samtools/1.20`, `/admin/software/samtools/samtools-1.20/bin` |
| ImageMagick | 7.1.2-33 | `/usersoftware/core001/common/RHEL_8/imagemagick/7.1.2-33/bin/magick`; personal module only, so `bin/loadTools.sh` prepends the path. |
| `/usr/bin/time` | absent | Peak RSS below was sampled with `ps` every 5 s. |
| vcf2maf | absent | Isoform override files vendored in `resources/vcf2maf/v1.6.12/`. |

## seqDNAcopy under R 4.5.1

The fork (soccin/seqDNAcopy 0.0.3) failed to install with

```
Error: object 'RangesList' is not exported by 'namespace:IRanges'
```

`RangesList` was removed from IRanges 2.x; `IRangesList` is the
replacement and is what `R/bams2counts.R` already calls. One line in
NAMESPACE; 0.0.4. Builds with g++ 12.2.0 from the R module's spack
toolchain.

## Measurements

Foreground on `isca006`, 2026-10-20, Proj_12257 pair
`s_C_71VKFN_N001_d01` (4.0 GB) x `s_C_TTH188_N001_d01` (5.0 GB), b37
exome, `BINSIZE=100`:

| Tool | Elapsed | Peak RSS |
|---|---|---|
| `getPairedCounts` | 176 s | 2.0 GB |
| `seqSegment` | 66 s | 2.7 GB |

## Validation runs

**Full run, 2026-10-20**, `../test01/`, control job 546394 on
`cmobic_cpu`; stage jobs 546419-546430 (FIX), 546510-546532 (COUNT),
548253-548272 (SEG), 561199 (POST), all on `cmobic_short`. 28 of 28
`COMPLETED`, 6 min 10 s wall. Per-stage Elapsed and MaxRSS in
`docs/SLURM_PORT.md`.

`MaxRSS` for the COUNT jobs ranged from 1.9 to 8.4 GB against
`--mem=8G` while the same tool peaks at 2.0 GB RSS in the foreground:
as PEMapper found, `jobacct_gather/cgroup` reports the cgroup
watermark including page cache, so a job streaming 9 GB of BAM fills
its cgroup with cache up to the cap. No job was OOM-killed.

**Failure tests, 2026-10-20**: `../test02/` (missing index, no job
submitted), `../test01/` (rerun refused), `../test03/` (control job
588179, one COUNT job cancelled by hand).
