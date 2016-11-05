args=commandArgs(trailing=T)
normalBam=args[1]
tumorBam=args[2]

tBase=gsub("indelRealigned_recal_","",gsub("___MD","",gsub(".bam","",basename(tumorBam))))
nBase=gsub("indelRealigned_recal_","",gsub("___MD","",gsub(".bam","",basename(normalBam))))

print(tBase)
print(nBase)

if(tBase==nBase){
    stop("tumor==normal")
}

sampleId=cc(tBase,"_",nBase)
print(sampleId)

library(Rsamtools)
library(DNAcopy)
library(pctGCdata)

source("SeqDNACopy/seqDNAcopy.R")

#There are 2 normals:
# CL_NL8533-2P  (this one was a "not real" line -- it was called CL_WD8533-2P when we ran it on CGH)
# CL_L070711   (preadipocyte cell line)

BINSIZE=10000

bb=bams2counts(normalBam,tumorBam,X=T)
out=seqsegment(bb,sampleid=sampleId,binSize=BINSIZE)

save(bb,out,file=cc(sampleId,"Bin",BINSIZE,".Rdata"),compress=T)

