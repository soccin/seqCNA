## METHODS

### seqCNA: Copy Number Variant Calling from Sequencing Data

Copy number changes are identified by comparing the counts in the target sample with counts in a selected normal sample. To minize noise it is critical to pick an optimal normal sample for each target/tumor sample. We start from a set of curated normal samples that have been processed and sequenced in the same way as the target. We first compare the insert size distribution between the target and the set of normals and select 10 normals as first pass set to use. For each of these normals we first count the number of reads in 100bp bins for both the target and the normals. Bins within a blacklist of problematic areas are removed. A GC-bias correction is done by first computing a scaling factor as the ratio of the sum of counts in normal versus the target and then normalized log2 ratio of target to normal is computed. A GC-bias correction is done using a weighted loess regression of the logR to the pctgc for each bin. The weighting factor is the log or the sum of the tumor and normal counts. The fitted curve is used to remove the GC-bias by scaling the raw target counts.

The scaled target/tumor counts and raw normal counts are then segmented using the CBS algorithm [Olshen 2004]. First any bins that have less then 35 counts (minBinCounts) in the normal sample are discarded. Then for the remaining bins a scaled log ratio of target to normal is compute again this time though adding a pseudo count of 1 to handle log(0). A weighting factor of log2(NormalCounts+1-minBincounts) is used and the standard DNAcopy::segment method is used to compute the segmentation with sdundo for undo.splits and an undo.SD=2. The mean segment values are then clustered/merged by merging any segments whose mean value differ by less than a threshold of 0.04. This reduces the low level noise often seen between chromosome segments. Finally the segment cluster which has the smallest logR/sqrt(numProbes) is set to logR zero to set the diploid level. These diploid segments are used subsequently to compute significantly altered segments.

This procedure is repeated for each normal in the set of choosen above. Then the normal with the smallest value of the sum of logR squared is picked as the optimal normal and the output from this pair is used as the copy number profile for the given target sample. To determine which segments are significantly altered we used the method from [Cheng 2015] the standard deviation of all probe values in all of the segments assigned to the diploid cluster is computed. This standard deviation is used in a 1 sample t-test for each segment comparing the segment mean to this standard deviation. The p-value is computed and then corrected using the FDR method. Segments are scored as significant if the FDR<0.05 and the absolute value of the log2 ratio is greater then 1 (abs fold change > 2).


### References:

- Olshen, AB, et al; Biostatistics 2004 Oct; 5(4): 557-572

- Cheng, D, et al; J Mol Diagn. 2015 May; 17(3): 251–264.


### Github links:

- (https://github.com/soccin/seqCNA)
- (https://github.com/soccin/seqDNAcopy)


