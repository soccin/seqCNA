#
# doWGSCNA.R (version 3.0.1)
#

require(Rsamtools)
require(DNAcopy)
require(pctGCdata)
require(Cairo)
require(stringr)

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

source(file.path(getSDIR(),"SeqDNACopy/seqDNAcopy.R"))

cArgs=commandArgs(trailing=T)


################################################################
#
# This code will parse command line args in the form of
#    KEY=VAL
# and sets
#    args[[KEY]]=VAL
#

# Set defaults first

args=list()
parseArgs=str_match(cArgs,"(.*)=(.*)")
apply(parseArgs,1,function(x){args[[str_trim(x[2])]]<<-str_trim(x[3])})

################################################################

normalBam=args$NORMAL
tumorBam=args$TUMOR

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
binSize=100*floor(binSize/quantile(nCounts,.25))
cat("adjusted binSize =",binSize,"\n")


out=seqsegment(bb,sampleid=sampleId,binSize=binSize)
out$param=list()
out$param$binSize=binSize
out$param$sampleId=sampleId

offset=c(0,tapply(bb$pos,bb$chrom,max))
gPos=cumsum(offset)/1e6

chromoLabels=c(seq(22),"X")
stagger=.1*((seq(chromoLabels) %% 2)-.5)

#CairoPDF(file=cc(sampleId,"Bin",binSize,".pdf"),width=11,height=8,bg="white",pointsize=12)

png(file=cc(sampleId,"Bin",binSize,".png"),
        type="cairo",
        width=1150,height=800,pointsize=20)

plot(out,xmaploc=T,ylim=c(-3,3))
abline(v=gPos,lty=2,col=8)
abline(h=c(-1,1),lty=2,col="#DDDDDD",lwd=1)
text(gPos[-len(gPos)]+offset[-1]/2e6,-1+stagger,chromoLabels,cex=.71)
text(0.5,1.5,
    paste(
        "BinSize =",
        formatC(out$param$binSize,format="d",big.mark=","),
        "NumBins =",
        formatC(nrow(out$data),format="d",big.mark=",")),
    pos=4,cex=1.414)

dev.off()

save(bb,out,file=cc(sampleId,"Bin",binSize,".Rdata"),compress=T)
