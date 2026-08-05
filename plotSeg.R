#
# plotSeg.R
#
# Standalone 1-panel plot of a seqSegment segmentation (_seqSeg.rda)
#
# usage: Rscript plotSeg.R [-b|--base] [-r|--redact] SEG.rda [OUT.png] [YLIM]
#
#    -b|--base     write the png to the current directory
#    -r|--redact   strip "__s_.*$" from the sample id in the title
#

cArgs=commandArgs(trailing=T)

BASE=any(cArgs %in% c("-b","--base"))
REDACT=any(cArgs %in% c("-r","--redact"))
cArgs=cArgs[!cArgs %in% c("-b","--base","-r","--redact")]

if(length(cArgs)<1) {
    cat("\n   usage: Rscript plotSeg.R [-b|--base] [-r|--redact] SEG.rda [OUT.png] [YLIM]\n\n")
    quit()
}

segFile=cArgs[1]
pFile=ifelse(length(cArgs)>1,cArgs[2],sub("\\.rda$",".png",segFile))
YLIM=ifelse(length(cArgs)>2,as.numeric(cArgs[3]),5)

if(BASE) pFile=basename(pFile)

suppressPackageStartupMessages(library(DNAcopy))

load(segFile) # loads: out

#
# The title is colnames(out$data)[3]; out$output$ID must match it
# or plot.DNAcopy will not draw the segment means
#

if(REDACT) {
    sampleId=sub("__s_.*$","",colnames(out$data)[3])
    colnames(out$data)[3]=sampleId
    out$output$ID=sampleId
}

#
# MAD of the log-ratios about their segment means
#

probe.seg.values=rep(out$output$seg.mean,out$output$num.mark)
global.mad=mad(out$data[,3]-probe.seg.values)

plot1Panels <- function(out) {

    plot(out,xmaploc=T,ylim=YLIM*c(-1,1),pt.cols=c("#B5D7E4","#BEBEBE"))

    abline(h=c(-1,1,log2(1.5)),lty=2,col="#333333",lwd=1)
    abline(h=global.mad*c(-1,1),lty=3,col=1)

    chrom.max=tapply(out$data$maploc,out$data$chrom,max)
    abline(v=c(0,cumsum(chrom.max)),lty=2,lwd=2,col="#666666")
    text(cumsum(chrom.max)-chrom.max/2,-YLIM,unique(out$data$chrom),cex=0.8)

    box()

}

#
# same 1150x800 layout as seqSegment, rendered at DPI instead of 72
#

DPI=150

png(file=pFile,type="cairo",width=1150/72,height=800/72,units="in",res=DPI,pointsize=16)
plot1Panels(out)
dev.off()

cat("#",pFile,"\n")
