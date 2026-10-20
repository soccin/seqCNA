# seqCNA 6.0.0: release checklist

Everything that must be done before `devs` is merged to `master` and
tagged 6.0.0. The port itself is merged into `devs` (bae5001). The
open issues are in `00.ISSUES.md`; the port spec and its decisions
(D1-D16) are in `attic/docs/00.SPEC_IRIS_PORT.md`.

## 1. R environment

- [ ] Install DNAcopy and pctGCdata in the personal R 4.5 library
      (`BiocManager::install("DNAcopy")`,
      `install.packages("pctGCdata")`, or
      `remotes::install_github("mskcc/pctGCdata")`). Spec D2.
- [ ] Delete the temporary copies `Rlib/DNAcopy` and `Rlib/pctGCdata`.
      They were copied from the trial library so the test jobs could
      run, and they hide the check in `00.SETUP.sh`.
- [ ] Rerun `00.SETUP.sh`. The switch to `devs` and the merge removed
      the built `Rlib/seqDNAcopy` (the old binary was tracked on
      `devs`), so the pipeline will not start until it is rebuilt.
      Check `Rlib/VERSIONS.txt` shows the packages coming from the
      personal library and seqDNAcopy 0.0.4.

## 2. Validation

- [ ] Rerun the full Proj_12257 test (as `../test01`) from a clean
      `devs` checkout with the personal library only. The 2026-10-20
      runs used a dirty tree and the temporary Rlib copies.
- [ ] Failure test not yet run: `scancel` the control job itself while
      stage jobs are queued (spec section 4, step 3). Expect `FAILED`
      in `00.RUNSTATUS.txt` and no stray jobs.
- [ ] Compare with a JUNO run of a previous seqCNA project
      (`00.ISSUES.md` #1, spec D12 and section 4 step 5). Diff
      `_seqSeg.seg`, `_seqSeg.out`, `bestMatches____out`,
      `GeneMatrix.csv`. Differences must be explained by package
      versions and stay below the calling thresholds.
- [ ] Run a mouse project (`mm10` relabel, `Mouse-Twist`) and an sWGS
      project (`BINSIZE=auto`), or keep them under "Not tested" in the
      CHANGELOG (`00.ISSUES.md` #2).
- [ ] Revisit the COUNT and SEG memory request (`-M 8G`) after the sWGS
      run (`00.ISSUES.md` #7).
- [ ] Record the new runs (job ids, sacct Elapsed and MaxRSS) in
      `docs/SLURM_PORT.md` and `docs/IRIS_CLUSTER_INFO.md`.

## 3. Docs

- [ ] Decision D15 in `attic/docs/00.SPEC_IRIS_PORT.md` says "so no
      tag"; the annotated tag `v5-juno` now exists on 3e6c18d. Fix
      D15, and name the tag next to the `juno` branch in `CHANGELOG.md`
      and `README.md`.
- [ ] Pooled normals (`IDT_Exome_v2_Normal`): add the placeholder entry
      in `docs/IRIS_CLUSTER_INFO.md` that D10 calls for; fill in the
      IRIS path in it and in `README.md` once the BAMs are copied
      (`00.ISSUES.md` #3).
- [ ] CHANGELOG: record the source and md5 of the two vendored isoform
      files (D8). They are in `resources/vcf2maf/README.md` only.
- [ ] `METHODS.md`: add software versions (R 4.5.1, DNAcopy 1.84.0,
      seqDNAcopy 0.0.4) if the method text should state them (spec
      section 3.6). It currently names no versions.
- [ ] Reread `CLAUDE.md`, `README.md`, `docs/*.doc` against the final
      code.

## 4. Release

- [ ] Bump `VERSION` from `5.dev` to `6.0.0` in `getPairedCounts.R`,
      `seqSegment.R`, `getGeneCalls.R` and `Version:` in `DESCRIPTION`
      (D15).
- [ ] CHANGELOG: `[Unreleased]` becomes `[6.0.0] - <date>`; move items
      under "Not tested" as validation closes them.
- [ ] Merge `devs` -> `master`; tag the release (existing release tags
      use `ver_X.Y.Z`; choose `ver_6.0.0` or `v6.0.0`). `PIPE.sh` and
      `getGITTag` use `git describe --tags`, so the tag is what the run
      records show.
- [ ] Push `devs`, `master`, the release tag and `v5-juno` to origin.
      Delete `feat/iris` once merged.

## 5. Cleanup

- [ ] Delete scratch: `/localscratch/soccin/seqcna-portref`
      (reference port clones, trial `Rlib45`),
      `/localscratch/soccin/seqcna-handrun`.
- [ ] Delete or keep the test directories `../test01`, `../test02`,
      `../test03` (the job ids in the docs refer to them).
- [ ] Remove this file at release, or move it to `attic/docs/`.

## After 6.0.0

Not required for the release; tracked in `00.ISSUES.md`:
`delivery.sh` for IRIS (#4), one genome table (#5), the unused
`samples` file (#6), status file on usage errors (#8). Superseded
tidyverse calls (`top_n`, `summarize_at`, `spread`) and the undefined
`first_UP()` in `convert_mouse_to_human.R` (spec section 2.3).
