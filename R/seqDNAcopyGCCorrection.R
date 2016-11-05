library(DNAcopy)

source("counts2logR.R")

doCNA<-function(zz,sampleid,minBinCount) {

    # weights
    # zwts <- log2(zz[,"rCountN"]+1-minBinCount)
    zwts <- log2(1/{1/(zz[,"rCountT"]+1) + 1/(zz[,"rCountN"]+1)})
    ii=which(zwts>0)

    zwts=zwts[ii]

    zchr <- zz[ii,"chrom"]
    zpos <- zz[ii,"maploc"]
    # log-ratio
    zlr <- zz[ii,"cnlr"]



    zcna <- CNA(zlr, zchr, zpos, sampleid=sampleid)
    zout <- segment(zcna, weights=zwts, undo.splits="sdundo",undo.SD=2)
    #zout <- segment(zcna)
    zout
}


# From BEDTOOLS multicov

COUNTFILE="counts.txt"
BAMFILES="bamList2"
#normal="Mtb_H37Rv_North_CN2015"
normal="s_Mtub_WT_neu_S1"


BAMS=scan(BAMFILES,"")

samples=gsub(".bam","",gsub(".*recal_","",BAMS))
counts=read.delim(COUNTFILE,header=F)

colnames(counts)=c("chrom","start","stop",samples)
counts$pos=(counts$start+counts$stop)/2e6

minBinCount=15


for(si in samples){
    if(si != normal){
        print(si)

        mat=data.frame(chrom=as.numeric(factor(counts$chrom)),
                        maploc=floor((counts$start+counts$stop)/2),
                        rCountN=counts[,normal],
                        rCountT=counts[,si])
        mat$keep=rep(1,nrow(mat))
        mat$keep[mat$rCountN<minBinCount]=0

        out=counts2logR(mat)


        zout=doCNA(out,si,minBinCount)
        plot(zout,xmaploc=T,ylim=c(-10,10))
        dev.copy2pdf(file=cc("dnaCopy",si,"_vs_",normal,".pdf"),width=11,height=8)
        write.xls(zout$output,file=cc("dnaCopy",si,"_vs_",normal,".seg"))
    }
}

