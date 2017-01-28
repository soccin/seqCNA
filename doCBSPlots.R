library(Cairo)
library(DNAcopy)

offset=NULL

rdat=dir(pattern=".Rdata")

Ns=len(rdat)
stats=data.frame(
    Tsample=rep("",Ns),
    Nsample=rep("",Ns),
    BinSize=rep(NA,Ns),
    NumBins=rep(NA,Ns)
    )

for(ii in seq(rdat)) {

    rfile=rdat[ii]
    print(rfile)
    load(rfile)

    if(is.null(offset)){
        offset=c(0,tapply(bb$pos,bb$chrom,max))
        gPos=cumsum(offset)/1e6
        chromoLabels=c(seq(22),"X")
        stagger=.1*((seq(chromoLabels) %% 2)-.5)
    }

    #png(file=cc(out$output$ID[1],".png"),
     #   type="cairo",
     #   width=2400, height=1600, bg="white", dpi=300)

    png(file=cc(out$output$ID[1],".png"),
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

    sInfo=strsplit(rfile,"___")[[1]]
    stats[ii,]=list(
        sInfo[1],
        gsub("_Bin.*","",sInfo[2]),
        as.numeric(out$param$binSize),
        nrow(out$data))
}

library(xlsx)
write.xlsx2(stats,file=cc("sWG_CNA__Stats__07260",DATE(),".xlsx"),row.names=F)
