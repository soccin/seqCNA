The results for the DNA copy number pipeline are now ready. You can
access them on the server at:

    https://bicdelivery.mskcc.org/project/${PROJECTNO}/r_001

in the folder

    seqcna

You will find the following files:

- *.pdf is a multipage PDF with the copy number profiles for each
  sample

- .seg file has the output of the DNAcopy CBS segmentation algorithm.
  This file can be loaded into IGV to view a genome-wide figure.

- GeneMatrix.csv, SegmentMatrix.csv, Proj_${PROJECTNO}___GeneTable.xlsx
  These files have gene-level and segment-level _significant_ calls,
  where segments are scored as significant if the FDR < 0.05 and
  abs(log2R) > 1

  * Proj_${PROJECTNO}___GeneTable.xlsx has significant gene calls with
    the spanning segments along with their FDR value and log2 ratios

  * SegmentMatrix.csv, GeneMatrix.csv have the gene-level significant
    calls in matrix form with log2 ratio values. SegmentMatrix.csv
    includes the data for the spanning segment. GeneMatrix.csv just has
    the gene info

Finally, if you would like the output BAM files for further analysis,
you can get them in the alignment folder.

If you have any questions let me know.

Nicholas Socci
Bioinformatics Core
MSKCC
soccin@mskcc.org
