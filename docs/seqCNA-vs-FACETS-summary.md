# seqCNA vs FACETS: Summary

**seqCNA** computes log2 tumor/normal coverage ratios, segments them with CBS, and reports
gene-level copy number calls. It is simple, robust, and handles cases FACETS cannot:
unmatched or pooled normals, mouse genomes, and mixed capture platforms.

**FACETS** jointly models the coverage log-ratio and B-allele frequencies at germline SNPs
to produce integer allele-specific copy number, tumor purity, ploidy, and LOH calls. It is
more informative but has more failure modes.

| | seqCNA | FACETS |
|---|---|---|
| Output | Continuous log2 ratio | Integer copy number |
| Purity / ploidy | Not estimated | Estimated |
| LOH detection | No | Yes |
| Matched normal required | No | Yes |
| Mouse genome support | Yes | No |
| Failure rate | Low | Non-trivial |

**Use seqCNA** when there is no matched normal, for mouse tumors, shallow WGS, or
cross-platform cohorts where a single consistent method is needed.

**Use FACETS** when you have a matched tumor/normal pair, need integer copy number, LOH,
or purity/ploidy estimates, or are integrating with TEMPO.

**Key seqCNA caveat:** Output is relative only -- no absolute copy number, no purity
correction, no LOH. A segment at log2R = -1 could represent one lost copy in a pure
diploid tumor or any number of other scenarios; the method cannot resolve this ambiguity.
It is explicitly a fallback for situations where FACETS cannot run.

**Key FACETS caveat:** Fails or degrades on low-purity samples, inbred mouse models (too
few heterozygous SNPs), mismatched bait sets, and highly aneuploid or genome-doubled
tumors where the diploid baseline is unreliable. Requires expert review; automated failure
flags do not catch all bad solutions.
