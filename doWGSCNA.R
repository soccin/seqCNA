#
# doWGSCNA.R (version 2.2)
#

library(Rsamtools)
library(DNAcopy)
library(pctGCdata)
library(Cairo)

source(file.path(getSDIR(),"SeqDNACopy/seqDNAcopy.R"))

################################################################
fixSampleNames<-function(x) {

    x=gsub("indelRealigned_recal_","",x)
    x=gsub("___MD","",x)
    x=gsub("^s_","",x)
    x=gsub(".bam$","",x)
    return(x)

}

getSDIR <- function(){
    args=commandArgs(trailing=F)
    TAG="--file="
    path_idx=grep(TAG,args)
    SDIR=dirname(substr(args[path_idx],nchar(TAG)+1,nchar(args[path_idx])))
    if(length(SDIR)==0) {
        return(getwd())
    } else {
        return(SDIR)
    }
}
################################################################

args=commandArgs(trailing=T)
normalBam=args[1]
tumorBam=args[2]

tBase=fixSampleNames(basename(tumorBam))
nBase=fixSampleNames(basename(normalBam))

print(tBase)
print(nBase)

if(tBase==nBase){
    stop("tumor==normal")
}

sampleId=cc(tBase,"_",nBase)
print(sampleId)


bb=bams2counts(normalBam,tumorBam,X=T)


binSize=50000
bins=bb$chrom + floor(bb$pos/binSize)*binSize/2^28
nCounts=tapply(bb$normal,bins,sum)
quantile(nCounts)
binSize=floor(binSize*100/quantile(nCounts,.25))
cat("adjusted binSize =",binSize,"\n")


out=seqsegment(bb,sampleid=sampleId,binSize=binSize)


offset=c(0,tapply(bb$pos,bb$chrom,max))
gPos=cumsum(offset)/1e6

chromoLabels=c(seq(22),"X")
stagger=.1*((seq(chromoLabels) %% 2)-.5)

CairoPDF(file=cc(sampleId,"Bin",binSize,".pdf"),width=11,height=8,bg="white",pointsize=12)

plot(out,xmaploc=T,ylim=c(-2,2))
abline(v=gPos,lty=2,col=8)
abline(h=c(-1,1),lty=3,col="#DDDDDD",lwd=2)
text(gPos[-len(gPos)]+offset[-1]/2e6,-1+stagger,chromoLabels,cex=.71)

dev.off()

save(bb,out,file=cc(sampleId,"Bin",binSize,".Rdata"),compress=T)
