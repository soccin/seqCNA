#
# doWGSCNA.R
#

VERSION="5.dev"

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


cArgs=commandArgs(trailing=T)
suppressPackageStartupMessages(require(stringr))

################################################################
#
# This code will parse command line args in the form of
#    KEY=VAL
# and sets
#    args[[KEY]]=VAL
#

# Set defaults first
#   require set to NULL

args=list(
    BINSIZE="auto",
    TUMOR=NULL,
    NORMAL=NULL,
    MAPLOC="false",
    GENOME="hg19",
    ODIR="."
    )

parseArgs=str_match(cArgs,"(.*)=(.*)")
dummy=apply(parseArgs,1,function(x){args[[str_trim(x[2])]]<<-str_trim(x[3])})

cat("\n\n###############################################################\n")
cat("#\n# Version:",VERSION,"\n")

if(any(sapply(args,is.null))) {
    cat("\n\nusage: doWGSCNA.R TUMOR=/path/tumor.bam NORMAL=/path/normal.bam BINSIZE=[auto] MAPLOC=[false] GENOME=[hg19]\n\n")
    missing=which(sapply(args,is.null))
    cat("missing required arg(s)\n\n   ")
    for(ii in missing){
        cat(names(args)[[ii]],"")
    }
    cat("\n\n")
    quit()
}

if(args$ODIR!="." & !dir.exists(args$ODIR)) {
    dir.create(args$ODIR,recursive=T)
}

################################################################
# Load rest of libraries after args checkout
#

suppressPackageStartupMessages(library(seqDNAcopy,lib.loc=file.path(getSDIR(),"Rlib")))
cat("# Version(seqDNAcopy):",sessionInfo()$otherPkgs$seqDNAcopy$Version,"\n")
#require(Rsamtools)
#require(DNAcopy)
suppressPackageStartupMessages(require(Cairo))
cat("#\n\n")

keys=sort(names(args))
for(key in keys) {
    cat("#",key,"=",args[[key]],"\n")
}

################################################################

normalBam=args$NORMAL
tumorBam=args$TUMOR

tBase=fixSampleNames(basename(tumorBam))
nBase=fixSampleNames(basename(normalBam))

if(tBase==nBase){
    stop("tumor==normal")
}

sampleId=cc(tBase,"_",nBase)
cat("# sampleId =",sampleId,"\n")

bb=bams2counts(normalBam,tumorBam,X=T,gbuild=args$GENOME)

if(args$BINSIZE=="auto") {
    binSize=50000
    bins=bb$chrom + floor(bb$pos/binSize)*binSize/2^28
    nCounts=tapply(bb$normal,bins,sum)
    quantile(nCounts)
    binSize=100*floor(binSize/quantile(nCounts,.25))
    cat("# adjusted binSize =",binSize,"\n")
} else {
    binSize=as.numeric(args$BINSIZE)
    cat("# manually set binsize =",binSize,"\n")
}

out=seqsegment(bb,sampleid=sampleId,binSize=binSize)
out$param=list()
out$param$binSize=binSize
out$param$sampleId=sampleId

offset=c(0,tapply(bb$pos,bb$chrom,max))
gPos=cumsum(offset)/1e6

chromoLabels=c(seq(22),"X")
stagger=.1*((seq(chromoLabels) %% 2)-.5)

png(file=file.path(args$ODIR,cc(sampleId,"Bin",binSize,".png")),
        type="cairo",
        width=1150,height=800,pointsize=20)

YLIM=3
if(args$MAPLOC=="true" & args$GENOME=="hg19") {

    plot(out,xmaploc=T,ylim=YLIM*c(-1,1))
    abline(v=gPos,lty=2,col=8)
    text(gPos[-len(gPos)]+offset[-1]/2e6,-1+stagger,chromoLabels,cex=.71)

} else {

    plot(out,xmaploc=F,ylim=YLIM*c(-1,1))

}

abline(h=c(-1,1),lty=2,col="#DDDDDD",lwd=1)
text(0.5,3.5,
    paste(
        "BinSize =",
        formatC(out$param$binSize,format="d",big.mark=","),
        "NumBins =",
        formatC(nrow(out$data),format="d",big.mark=",")),
    pos=4,cex=1.414)

dev.off()

save(bb,out,file=file.path(args$ODIR,cc(sampleId,"Bin",binSize,".Rdata")),compress=T)
