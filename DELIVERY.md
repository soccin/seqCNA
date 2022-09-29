In the folder

    seqcna

You will find the following files:

- *.pdf which is a multipage PDF which the copy number profiles for
each sample

- .seg file which has the output of the DNAcopy CBS segmentation
algorithm. This file can be loaded for view in IGV
a genome wide figure from IGV

- GeneMatrix.csv, SegmentMatrix.csv, ${PROJECTNO}___GeneTable.xlsx
This files have gene level and segment level _significant_ calls where
segments scored to be significant as having an FDR < 0.05 and abs(log2R)>1

  * ${PROJECTNO}___GeneTable.xlsx has significant gene calls with the
    spanning segments along with their FDR value and log2 Ratios

  * SegmentMatrix.csv, GeneMatrix.csv have the gene level significant
    in a matrix form with log2 Ratio values. SegmentMatrix.csv includes
    the data for the spanning segment. GeneMatrix.csv just have the gene
    info

Finally if you would like to output the output BAM files for further
analysis you can get them in the alignment folder.

If you have any questions let me know.

Nicholas Socci
Bioinformatics Core
MSKCC
soccin@mskcc.org
