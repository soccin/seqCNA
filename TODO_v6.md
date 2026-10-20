# seqCNA 6.0.0: release checklist

Open items before `rel/v6.0.0` is merged to `master` and tagged 6.0.0.
Open issues, including the ones not required for the release, are in
`00.ISSUES.md`.

## 1. Validation

- [ ] Compare with a JUNO run of a previous seqCNA project
      (`00.ISSUES.md` #1, spec D12 and section 4 step 5). Diff
      `_seqSeg.seg`, `_seqSeg.out`, `bestMatches____out`,
      `GeneMatrix.csv`. Differences must be explained by package
      versions and stay below the calling thresholds. In progress in
      `../test20`; record its job ids and sacct Elapsed and MaxRSS in
      `docs/SLURM_PORT.md` and `docs/IRIS_CLUSTER_INFO.md`.

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
