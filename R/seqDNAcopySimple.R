library(DNAcopy)

doCNA<-function(counts,normal,tumor,minBinCount=15) {

    zz = counts[counts[,normal]>minBinCount,c("chrom","pos",normal,tumor)]
    colnames(zz)[3]="normal"
    colnames(zz)[4]="tumor"

    tscl <- sum(zz[,"normal"])/sum(zz[,"tumor"])
    zchr <- zz[,"chrom"]
    zpos <- zz[,"pos"]
    # log-ratio
    zlr <- log2(zz[,"tumor"]*tscl+1) - log2(zz[,"normal"]+1)
    # weights
    # zwts <- log2(1/{1/(zz[,"tumor"]+1) + 1/(zz[,"normal"]+1)})
    zwts <- log2(zz[,"normal"]+1-minBinCount)

    zcna <- CNA(zlr, zchr, zpos, sampleid=tumor)
    zout <- segment(zcna, weights=zwts)
    zout
}


# From BEDTOOLS multicov

COUNTFILE="counts.txt"
BAMFILES="bamList"
BAMS=scan(BAMFILES,"")

samples=gsub(".bam","",gsub(".*recal_","",BAMS))
counts=read.delim(COUNTFILE,header=F)

colnames(counts)=c("chrom","start","stop",samples)
counts$pos=(counts$start+counts$stop)/2e6

minBinCount=15

normal="Mtb_H37Rv_North_CN2015"

for(si in samples){
    if(si != normal){
        print(si)
        zout=doCNA(counts,normal,si)
        plot(zout,xmaploc=T)
        dev.copy2pdf(file=cc("dnaCopy",si,"_vs_",normal,".pdf"),width=11,height=8)
    }
}

