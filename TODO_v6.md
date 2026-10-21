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
- [x] Pin versions under R 4.5.1: a lock file of the installed
      packages and their dependencies, installed into `Rlib/` by
      `00.SETUP.sh`; one launcher restricts R to `Rlib` and stops on a
      version mismatch; samtools and ImageMagick versions recorded.
- [x] Small functional test, run twice (identical output with the
      seed): test31, test32 (`docs/SLURM_PORT.md`).
- [ ] The Proj_12257 rerun in `../../CNA` (user).

## 2. Release

- [ ] Merge `feat/delivery` (worktree `../seqCNA-delivery`) into
      `rel/v6.0.0` (`00.ISSUES.md` #4).
- [ ] Bump `VERSION` from `5.dev` to `6.0.0` in `getPairedCounts.R`,
      `seqSegment.R`, `getGeneCalls.R` and `Version:` in `DESCRIPTION`
      (D15).
- [ ] CHANGELOG: add the release date to the `[6.0.0]` heading; move
      items under "Not tested" as validation closes them.
- [ ] Merge `rel/v6.0.0` -> `master`; tag the release (existing
      release tags use `ver_X.Y.Z`; choose `ver_6.0.0` or `v6.0.0`).
      `PIPE.sh` and `getGITTag` use `git describe --tags`, so the tag
      is what the run records show.
- [ ] Push `devs`, `master`, the release tag and `v5-juno` to origin.
      Delete `feat/iris` once merged.
- [ ] Remove this file, or move it to `attic/docs/`.
