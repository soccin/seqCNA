counts2logR <-function(mat,f=0.2){

    out <- mat[mat$keep==1,]
    out$gcpct <- rep(NA_real_, nrow(out))
    # load gc data
    load("H37RvCOgcpct.rda")
    # loop thru chromosomes
    nchr <- max(mat$chrom) # IMPACT doesn't have X so only 22
    for (i in 1:nchr) {
        ii <- out$chrom==i
        jj <- ceiling((out$maploc[ii]-450)/100)
        jj[jj < 1] <- 1 # fix issues (maploc < 500) with chr17
        out$gcpct[ii] <- gcpctdb[[i]][jj]
    }
    ##### log-ratio with gc correction and maf log-odds ratio steps
    chrom <- out$chrom
    maploc <- out$maploc
    rCountN <- out$rCountN
    rCountT <- out$rCountT
    gcpct <- out$gcpct

    # compute gc bias
    ncount <- tapply(rCountN, gcpct, sum)
    tcount <- tapply(rCountT, gcpct, sum)
    pctgc <- as.numeric(names(ncount))
    tscl <- sum(ncount)/sum(tcount)
    gcb <- lowess(pctgc, log2(tcount*tscl)-log2(ncount), f=f)
    jj <- match(gcpct, gcb$x)
    gcbias <- gcb$y[jj]

    # compute cn log-ratio (gc corrected)
    cnlr <- log2(1+rCountT*tscl) - log2(1+rCountN) - gcbias

    out$cnlr <- out$gcbias <- rep(NA_real_, nrow(out))
    out$gcbias <- gcbias
    out$cnlr <- cnlr

    out

}
