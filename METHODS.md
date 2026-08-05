# METHODS

## seqCNA: Copy Number Variant Calling from Sequencing Data

Copy number changes are identified by comparing the counts in the target sample with counts in a selected normal sample. To minimize noise it is critical to pick an optimal normal sample for each target/tumor sample. We start from a set of curated normal samples that have been processed and sequenced in the same way as the target. We first compare the insert size distribution between the target and the set of normals and select 10 normals as a first pass set to use. For each of these normals we first count the number of reads in 100bp bins for both the target and the normals. Bins within a blacklist of problematic areas are removed. The GC-bias correction starts by computing a scaling factor as the ratio of the sum of counts in the normal versus the target, and then a normalized log2 ratio of target to normal is computed for each GC percentage level. The GC-bias is estimated using a weighted loess regression of the logR to the pctgc. The weighting factor is the log2 of the sum of the tumor and normal counts. The fitted curve is used to remove the GC-bias by scaling the raw target counts.

The scaled target/tumor counts and raw normal counts are then segmented using the CBS algorithm [Olshen 2004]. First, any bins that have 35 or fewer counts (minBinCount) in the normal sample are discarded. Then for the remaining bins a scaled log ratio of target to normal is computed again, this time adding a pseudo count of 1 to handle log(0). A weighting factor of log2(NormalCounts+1-minBinCount) is used and the standard DNAcopy::segment method is used to compute the segmentation with undo.splits="sdundo" and undo.SD=2. The mean segment values are then clustered/merged by merging any segments whose mean values differ by less than a threshold of 0.04. This reduces the low-level noise often seen between chromosome segments. Finally, the segment cluster which has the smallest abs(logR)/sqrt(numProbes) is set to logR zero to set the diploid level. These diploid segments are used subsequently to compute significantly altered segments.

This procedure is repeated for each normal in the set chosen above. Then the normal with the smallest RMS logR, sqrt(mean(logR^2)), is picked as the optimal normal and the output from this pair is used as the copy number profile for the given target sample. To determine which segments are significantly altered we used the method from [Cheng 2015]: the mean and standard deviation of all probe values in all of the segments assigned to the diploid cluster are computed. For each segment a two-sided normal test is then used to compare the segment mean to the diploid cluster mean, using the diploid standard deviation as the scale. The resulting p-values are corrected for multiple testing using the FDR (Benjamini-Hochberg) method. Segments are scored as significant if the FDR < 0.05 and the absolute value of the log2 ratio is greater than 1 (absolute fold change > 2).


## References:

- Olshen, AB, et al; Biostatistics 2004 Oct; 5(4): 557-572

- Cheng, D, et al; J Mol Diagn. 2015 May; 17(3): 251–264.


## Links to source code

- https://github.com/soccin/seqCNA
- https://github.com/soccin/seqDNAcopy
