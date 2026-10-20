# seqCNA 6.0.0: release checklist

Open items before `rel/v6.0.0` is merged to `master` and tagged 6.0.0.
`rel/v6.0.0` is `master` with `devs` merged in; the docs pass is on
`docs/rel-v6` off it. Open issues are in `00.ISSUES.md`; the port spec
and its decisions (D1-D16) are in `attic/docs/00.SPEC_IRIS_PORT.md`.

## 1. Validation

- [ ] Failure test not yet run: `scancel` the control job itself while
      stage jobs are queued (spec section 4, step 3). Expect `FAILED`
      in `00.RUNSTATUS.txt` and no stray jobs.
- [ ] Compare with a JUNO run of a previous seqCNA project
      (`00.ISSUES.md` #1, spec D12 and section 4 step 5). Diff
      `_seqSeg.seg`, `_seqSeg.out`, `bestMatches____out`,
      `GeneMatrix.csv`. Differences must be explained by package
      versions and stay below the calling thresholds. Running in
      `../test20`; record its job ids and sacct Elapsed and MaxRSS in
      `docs/SLURM_PORT.md` and `docs/IRIS_CLUSTER_INFO.md`.
- [ ] Run a mouse project (`mm10` relabel, `Mouse-Twist`) and an sWGS
      project (`BINSIZE=auto`), or keep them under "Not tested" in the
      CHANGELOG (`00.ISSUES.md` #2).
- [ ] Revisit the COUNT and SEG memory request (`-M 8G`) after the sWGS
      run (`00.ISSUES.md` #7).

## 2. Release

- [ ] Merge `docs/rel-v6` into `rel/v6.0.0`.
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

## 3. Cleanup

- [ ] Delete scratch: `/localscratch/soccin/seqcna-portref`
      (reference port clones, trial `Rlib45`),
      `/localscratch/soccin/seqcna-handrun`,
      `/localscratch/soccin/seqcna-bug03`.
- [ ] Delete or keep the test directories `../test01` to `../test03`
      and `../test10` to `../test14` (the job ids in the docs refer to
      them).
- [ ] Remove this file at release, or move it to `attic/docs/`.

## After 6.0.0

Not required for the release; tracked in `00.ISSUES.md`:
`delivery.sh` for IRIS (#4), one genome table (#5), status file on
usage errors (#8), METHODS.md rewrite (#9), pooled normals on IRIS
(#3). Superseded tidyverse calls (`top_n`, `summarize_at`, `spread`)
and the undefined `first_UP()` in `convert_mouse_to_human.R` (spec
section 2.3).
