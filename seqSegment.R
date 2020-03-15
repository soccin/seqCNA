#
# getPairedCounts.R
#

VERSION="5.dev"

SNAME=Sys.getenv("SNAME")
SDIR=Sys.getenv("SDIR")
SDIR=ifelse(SDIR=="","seqCNA",SDIR)
source(file.path(SDIR,"include/parseArgs.R"))
source(file.path(SDIR,"include/tools.R"))
source(file.path(SDIR,"include/clusterSegs.R"))

args=list(
    COUNTS=NULL,
    ODIR=".",
    BINSIZE="auto",
    MINBINCOUNT=35,
    YLIM=5
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

library(RJSONIO)
suppressPackageStartupMessages(library(seqDNAcopy))
cat("# Version(seqDNAcopy):",packageDescription("seqDNAcopy")$Version,"\n")
keys=sort(names(args))
for(key in keys) {
    cat("#",key,"=",args[[key]],"\n")
}

################################################################

load(args$COUNTS)
sampleId=counts$param$sampleId
bb=counts$bb
cArgs=counts$args


if(tolower(args$BINSIZE)=="auto") {
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

undo.SD=2
out=seqsegment(bb,sampleid=sampleId,
                binSize=binSize,minBinCount=args$MINBINCOUNT,
                undo.splits="sdundo",undo.SD=undo.SD)

clusterThreshold=0.04
out.seqseg=out
out=clustersegs(out,threshold=clusterThreshold)

out$param=list()
out$param$binSize=binSize
out$param$sampleId=sampleId
out$param$countFile=args$COUNTS
out$param$arg.binsize=args$BINSIZE
out$param$minBinCount=args$MINBINCOUNT
out$param$undo=list(undo.splits="sdundo",undo.SD=undo.SD)
out$param$genome=cArgs$GENOME
out$param$clusterThreshold=clusterThreshold

#
# probe.seg.values is the segment means projected
# back on to probe space
#

probe.seg.values=double(nrow(out$dat))
probe.cluster.number=integer(nrow(out$dat))
for(ii in seq(nrow(out$output))) {
    probe.seg.values[out$segRows[ii,1]:out$segRows[ii,2]]=out$output$seg.mean[ii]
    probe.cluster.number[out$segRows[ii,1]:out$segRows[ii,2]]=out$output$cluster[ii]
}

global.mad=mad(out$dat[,3]-probe.seg.values)
rms.derivative.noise=sqrt(mean(diff(out$dat[,3])^2))
sum.logr.sq=sum(out$dat[,3]^2)
rms.logr=sqrt(mean(out$dat[,3]^2))
rms.logr.flat=sqrt(mean((out$dat[,3]-probe.seg.values)^2))
frac.logR.ltNeg2=mean(out$dat[,3] < -2)

#
# mean of diploid cluster and sd for dmp-style signifcance
#

cl0.mean <- mean(out$data[probe.cluster.number==out$cluster$diploidClusterNum,3])
cl0.sd <- sd(out$data[probe.cluster.number==out$cluster$diploidClusterNum,3])
pVal <- 2*pnorm(abs(out$output$seg.mean-cl0.mean),sd=cl0.sd,lower.tail=F)
FDR <- p.adjust(pVal,"fdr")

numSegments=nrow(out$output)

save(out,file=file.path(args$ODIR,paste0(sampleId,"_seqSeg",".rda")),compress=T)
output=out$output
output$loc.start=output$loc.start*1.e6
output$loc.end=output$loc.end*1.e6
output$pValue=pVal
output$FDR=FDR
write.table(output,
            file=file.path(args$ODIR,paste0(sampleId,"_seqSeg",".seg")),
            sep="\t",eol="\n",quote=F,row.names=F)


##################################################################################
##################################################################################

outFile=file.path(args$ODIR,paste0(sampleId,"_seqSeg",".out"))
writeVariable<-function(varName) {
    cat(varName,"=",get(varName),"\n",file=outFile,append=T)
}
cat("##############################################\n",file=outFile)

writeVariable("VERSION")
GITTag=getGITTag(SDIR)
writeVariable("GITTag")
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
writeVariable("undo.SD")
writeVariable("global.mad")
writeVariable("rms.derivative.noise")
writeVariable("sum.logr.sq")
writeVariable("rms.logr")
writeVariable("rms.logr.flat")
writeVariable("frac.logR.ltNeg2")

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

YLIM=args$YLIM
plot2Panels <- function(out) {
    par(mfrow=c(2,1))

    plot(out,xmaploc=F,ylim=YLIM*c(-1,1),pt.cols=c("#B5D7E4","#BEBEBE"))

    abline(h=c(-1,1,log2(1.5)),lty=2,col="#333333",lwd=1)
    abline(h=global.mad*c(-1,1),lty=3,col=1)

    text(0.5,YLIM+.5-1,
        paste(
            "NumBins =",
            formatC(nrow(out$data),format="d",big.mark=","),
            "MAD =",
            formatC(global.mad,format="f"),
            "nSegs =",formatC(numSegments,format="d",big.mark=","),
            "RMSD.Segs.noX =",formatC(RMSD.seg.mean.AUTO,format="f")
            ),
        pos=4,cex=1.12)

    text(0.5,YLIM+.5-1.4,
        paste(
            "sum.logr.sq =",
            formatC(sum.logr.sq,format="f"),
            "rms.logr =",
            formatC(rms.logr,format="f"),
            "rms.logr.flat =",
            formatC(rms.logr.flat,format="f")
            ),
        pos=4,cex=1.12)

    plot(out,xmaploc=F,pt.cols=c("#B5D7E4","#BEBEBE"))
    abline(h=c(-2,2),lty=2,col="#333333",lwd=1)
}

plot1Panels <- function(out) {

    plot(out,xmaploc=T,ylim=YLIM*c(-1,1),pt.cols=c("#B5D7E4","#BEBEBE"))

    abline(h=c(-1,1,log2(1.5)),lty=2,col="#333333",lwd=1)
    abline(h=global.mad*c(-1,1),lty=3,col=1)

    chrom.max=tapply(out$data$maploc,out$data$chrom,max)
    abline(v=c(0,cumsum(chrom.max)),lty=2,lwd=2,col="#666666")
    text(cumsum(chrom.max)-chrom.max/2,-YLIM,unique(out$data$chrom),cex=0.8)



    label1=paste(
            "NumBins =",
            formatC(nrow(out$data),format="d",big.mark=","),
            "MAD =",
            formatC(global.mad,format="f"),
            "nSegs =",formatC(numSegments,format="d",big.mark=","),
            "RMSD.Segs.noX =",formatC(RMSD.seg.mean.AUTO,format="f")
            )

    label2=paste(
            "sum.logr.sq =",
            formatC(sum.logr.sq,format="f"),
            "rms.logr =",
            formatC(rms.logr,format="f"),
            "rms.logr.flat =",
            formatC(rms.logr.flat,format="f")
            )

    maxLabelLen=max(strwidth(label1,cex=1.12),strwidth(label2,cex=1.12))
    THEIGHT=strheight("XXX",cex=1.12)

    rect(-10,YLIM+.5-1+THEIGHT,maxLabelLen+0.5,YLIM+.5-1.4-THEIGHT,col="white",border=NA)

    rect(0,YLIM+.5-1+.1,maxLabelLen+0.5,YLIM+.5-1.4-.1,col="white",border=NA)

    text(0.5,YLIM+.5-1,label1,pos=4,cex=1.12)

    text(0.5,YLIM+.5-1.4,label2,pos=4,cex=1.12)

    box()

}

pFile=file.path(args$ODIR,paste0(sampleId,"_seqSeg",".png"))

# png(file=pFile,type="cairo",height=1150,width=800,pointsize=16)
# plot2Panels(out)

png(file=pFile,type="cairo",width=1150,height=800,pointsize=16)
plot1Panels(out)

dev.off()

