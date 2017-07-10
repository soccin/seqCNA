#
# getPairedCounts.R
#

VERSION="5.dev"

SNAME=Sys.getenv("SNAME")
SDIR=Sys.getenv("SDIR")
SDIR=ifelse(SDIR=="",".",SDIR)
source(file.path(SDIR,"include/parseArgs.R"))
source(file.path(SDIR,"include/tools.R"))

args=list(
    COUNTS=NULL,
    ODIR=".",
    BINSIZE="auto"
    )

cat("###############################################################\n")
cat("#\n# Version:",VERSION,"\n")

args=parseArgs(args)

if(is.null(args)) {

    cat(readLines(file.path(SDIR,"docs",paste0(SNAME,".doc"))),sep="\n")
    cat("\n\n")
    quit()

}

if(args$ODIR!="." & !dir.exists(args$ODIR)) {
    dir.create(args$ODIR,recursive=T)
}

################################################################
# Load rest of libraries after args checkout
#

suppressPackageStartupMessages(library(seqDNAcopy))
cat("# Version(seqDNAcopy):",sessionInfo()$otherPkgs$seqDNAcopy$Version,"\n")
keys=sort(names(args))
for(key in keys) {
    cat("#",key,"=",args[[key]],"\n")
}

################################################################

load(args$COUNTS)
sampleId=dat[[1]]$sampleId
bb=dat[[1]]$bb
cArgs=dat[[1]]$args

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
out$param$countFile=args$COUNTS
out$param$arg.binsize=args$BINSIZE

#
# probe.seg.values is the segment means projected
# back on to probe space
#

probe.seg.values=double(nrow(out$dat))
for(ii in seq(nrow(out$output))) {
    probe.seg.values[out$segRows[ii,1]:out$segRows[ii,2]]=out$output$seg.mean[ii]
}
global.mad=mad(out$dat[,3]-probe.seg.values)
rms.derivative.noise=sqrt(mean(diff(out$dat[,3])^2))
numSegments=nrow(out$output)

png(file=file.path(args$ODIR,cc(sampleId,"seqSeg",".png")),
        type="cairo",
        width=1150,height=800,pointsize=20)

YLIM=3
plot(out,xmaploc=F,ylim=YLIM*c(-1,1),pt.cols=c("#B5D7E4","#BEBEBE"))

abline(h=c(-1,1,log2(1.5)),lty=2,col="#333333",lwd=1)
abline(h=global.mad*c(-1,1),lty=3,col=1)

text(0.5,YLIM+.5-1,
    paste(
        "BinSize =",
        formatC(out$param$binSize,format="d",big.mark=","),
        "NumBins =",
        formatC(nrow(out$data),format="d",big.mark=","),
        "Global.MAD =",
        formatC(global.mad,format="f"),
        "NumSegs =",formatC(numSegments,format="d",big.mark=",")
        ),
    pos=4,cex=1.2)

text(0.5,YLIM+.5-1.25,
    paste(
        "max|segMean| Auto =",
        formatC(max(abs(out$output$seg.mean)),format="f"),
        "RMSD(seg.mean) =",
        formatC(sqrt(mean(out$output$seg.mean^2)),format="f")
        ),
    pos=4,cex=1.2)

dev.off()

save(out,file=file.path(args$ODIR,cc(sampleId,"seqSeg",".rda")),compress=T)
output=out$output
output$loc.start=output$loc.start*1.e6
output$loc.end=output$loc.end*1.e6
write.table(output,
            file=file.path(args$ODIR,cc(sampleId,"seqSeg",".seg")),
            sep="\t",eol="\n",quote=F,row.names=F)


##################################################################################
##################################################################################

outFile=file.path(args$ODIR,cc(sampleId,"seqSeg",".out"))
writeVariable<-function(varName) {
    cat(varName,"=",get(varName),"\n",file=outFile,append=T)
}
cat("##############################################\n",file=outFile)

writeVariable("VERSION")
writeVariable("sampleId")
countFile=args$COUNTS
writeVariable("countFile")
genome=cArgs$GENOME
writeVariable("genome")
arg.binsize=args$BINSIZE
writeVariable("arg.binsize")
writeVariable("binSize")
numBins=nrow(out$data)
writeVariable("numBins")
writeVariable("global.mad")
writeVariable("rms.derivative.noise")
numSegments=nrow(output)
writeVariable("numSegments")

XChr=NA
if(cArgs$GENOME %in% c("hg19","b37"))
    XChr=23
if(cArgs$GENOME %in% c("mm10"))
    XChr=20

RMSD.seg.mean=sqrt(mean(output$seg.mean^2))
autoII=output$chrom!=XChr
numSegments.AUTO=length(which(autoII))

RMSD.seg.mean.AUTO=sqrt(mean(output$seg.mean[autoII]^2))

wiAuto=output$num.mark[autoII]/sum(output$num.mark[autoII])

wRMSD.seg.mean.AUTO=sqrt(sum( wiAuto*(output$seg.mean[autoII]^2) ))
max.abs.seg.mean.AUTO=max(abs(output$seg.mean[autoII]))

writeVariable("numSegments.AUTO")
writeVariable("RMSD.seg.mean")
writeVariable("RMSD.seg.mean.AUTO")
writeVariable("wRMSD.seg.mean.AUTO")
writeVariable("max.abs.seg.mean.AUTO")

xII=output$chrom==XChr
X.seg.mean.max=max(output$seg.mean[xII])
writeVariable("X.seg.mean.max")
X.seg.mean.min=min(output$seg.mean[xII])
writeVariable("X.seg.mean.min")

wiX=output$num.mark[xII]/sum(output$num.mark[xII])
X.seg.mean.avg=sum(wiX*output$seg.mean[xII])
writeVariable("X.seg.mean.avg")

##########################################
#
# Cluster segments
#

# segs=out$output
# nsegs=nrow(segs)

# cnlr=out$data[is.finite(out$data[,3]),3]

# segid <- rep(rank(segs$seg.mean, ties.method="random"), segs$num.mark)
# segid <- sort(segid)

# segs=segs[order(segs$seg.mean),]

# ocnclust<-1:nsegs
# ocnlevels<-segs$seg.mean

#     while ((length(ocnlevels) > 1) && (min(diff(ocnlevels)) < 0.04)) {
#         j <- which.min(diff(ocnlevels))
#         ocnlevels <- ocnlevels[-j]
#         ocnclust[ocnclust>j] <- ocnclust[ocnclust>j] - 1
#         segid[segid > j] <- segid[segid > j] - 1
#         ocnlevels[j] <- median(cnlr[segid==j])
#     }
