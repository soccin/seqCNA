# seqCNA 6.0.0: release checklist

Open items before `rel/v6.0.0` is merged to `master` and tagged 6.0.0.
Open issues, including the ones not required for the release, are in
`00.ISSUES.md`.

## 1. From the JUNO comparison

`00.ISSUES.md` #1 is closed (DNAcopy weights fix). Notes in
`../test20/NOTES_JUNO_IRIS_COMP.md`.

- [x] chrX: no chrX gene calls in any assay. `getGeneCalls.R` keeps
      the autosomes only, by genome, for segments and annotation. X
      stays in segmentation, `rms.logr`, the IGV `.seg` and the plots.
- [x] Seed: `seqSegment SEED=` with a fixed default, `set.seed()`
      before `seqsegment()`, written to the `.out` file.
- [x] Pin versions under R 4.5.1: a lock file of R and the packages
      that set the output (DNAcopy, pctGCdata, Rsamtools, data.table,
      seqDNAcopy); one launcher stops on a version mismatch; samtools
      and ImageMagick versions recorded.
- [x] Small functional test, run twice (identical output with the
      seed): test31, test32 (`docs/dev/SLURM_PORT.md`).
- [x] The Proj_12257 rerun in `../../CNA` (user), `COMPLETED`
      2026-10-21, at e19e87b.

## 2. Release

- [x] Merge `feat/delivery` into `devs` (ec6a5bc). `delivery.sh` is
      ported (`00.ISSUES.md` #4) and tested with a delivery of
      Proj_12257.
- [x] Merge the docs branches into `devs`: `docs/reorg`,
      `docs/methods-output` and `docs/cleanup-01`, stacked in that
      order.
- [x] Merge `devs` into `rel/v6.0.0` (a2f5cc7).
- [x] Bump `VERSION` from `5.dev` to `6.0.0` in `getPairedCounts.R`,
      `seqSegment.R`, `getGeneCalls.R` and `Version:` in `DESCRIPTION`
      (D15).
- [x] CHANGELOG: add the release date to the `[6.0.0]` heading; move
      items under "Not tested" as validation closes them.
- [ ] Merge `rel/v6.0.0` -> `master`; tag the release `v6.0.0`
      (decided 2026-10-21), the same name as `VERSION.md`: the delivery
      email links to `docs/output.md` and `docs/methods.md` at that
      tag. `PIPE.sh` and `getGITTag` use `git describe --tags`, so the
      tag is what the run records show.
- [ ] Push `devs`, `master`, the release tag and `v5-juno` to origin.
      Delete `feat/iris` once merged.
- [x] Remove this file, or move it to `attic/docs/`.
