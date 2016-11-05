# library(Rsamtools)
# library(RSQLite)
# library(DNAcopy)

chromRange <- function(i, chr="") {
  eval(parse(text=paste('RangesList("', chr, i,'"=IRanges(start=0, end=536870912L))',sep="")))
}

bam2fragments <- function(bamFile, X=FALSE) {
    # get position, mate position and insert size (fragment length)
    what <- c("pos","mpos","isize")
    # get first mate from properly paired reads which pass QC
    flag=scanBamFlag(isNotPassingQualityControls=FALSE, isPaired=TRUE, isFirstMateRead=TRUE, isDuplicate=FALSE)
    bam <- list()
    # autosomes
    for(i in 1:22) {
        which <- chromRange(i)
        param <- ScanBamParam(flag=flag, which = which, what = what)
        # scan the data
        bam[[i]] <- scanBam(bamFile, param=param)[[1]]
    }
    if (X) {
        # chromosome X
        which <- chromRange("X")
        param <- ScanBamParam(flag=flag, which = which, what = what)
        # scan the data
        bam[[23]] <- scanBam(bamFile, param=param)[[1]]
    }
    bam
}

# fragments to counts binned into intervals of size binSize
fragments2counts <- function(nfmid, tfmid, binSize=100) {
    x <- table(floor(tfmid/binSize))
    y <- table(floor(nfmid/binSize))
    out <- list()
    bins <-  union(names(x),names(y))
    out$pos <- (as.numeric(bins)*binSize + 0.5*binSize)/1e6
    out$normal <- out$tumor <- rep(0,length(bins))
    ii <- match(names(y), bins)
    out$normal[ii] <- y
    ii <- match(names(x), bins)
    out$tumor[ii] <- x
    out
}

# convert the fragment counts to a data matrix
fragments2datamatrix <- function(nbam, tbam, binSize=100, iSizeLim=c(76,500)) {
    # normal fragment midpoints
    nfmid <- lapply(nbam, function(x, ll, ul) {
        x$isize <- abs(x$isize)
        (pmin(x$pos, x$mpos) + x$isize/2)[x$isize >= ll & x$isize <= ul]
    }, iSizeLim[1], iSizeLim[2])
    # tumor fragment midpoints
    tfmid <- lapply(tbam, function(x, ll, ul) {
        x$isize <- abs(x$isize)
        (pmin(x$pos, x$mpos) + x$isize/2)[x$isize >= ll & x$isize <= ul]
    }, iSizeLim[1], iSizeLim[2])
    # bin the mid points into interval of size binSize
    binnedcounts <- list()
    nchr <- length(nbam)
    for(i in 1:nchr) binnedcounts[[i]] <- fragments2counts(nfmid[[i]], tfmid[[i]], binSize=binSize)
    # total number of bins with nonzero count
    nbins <- unlist(lapply(binnedcounts, function(x) {length(x$pos)}))
    # convert to matrix
    fcounts <- matrix(0, sum(nbins), 4)
    colnames(fcounts) <- c("chrom","pos","normal","tumor")
    fcounts[,1] <- rep(1:nchr, nbins)
    for(i in 1:nchr) {
        ii <- fcounts[,1]==i
        fcounts[ii,2] <- binnedcounts[[i]]$pos
        fcounts[ii,3] <- binnedcounts[[i]]$normal
        fcounts[ii,4] <- binnedcounts[[i]]$tumor
    }
    # return data
    fcounts
}

# segment the data
bam2segment <- function(nBamFile, tBamFile, sampleid, minBinCount=15, binSize=100, iSizeLim=c(76,500), X=FALSE) {
    nbam <- bam2fragments(nBamFile, X)
    tbam <- bam2fragments(tBamFile, X)
    fcounts <- fragments2datamatrix(nbam, tbam, binSize, iSizeLim)
    
    # counts to segments
    zz <- fcounts[fcounts[,3]>minBinCount,]
    # scale the read counts
    tscl <- sum(zz[,"normal"])/sum(zz[,"tumor"])
    zchr <- zz[,"chrom"]
    zpos <- zz[,"pos"]
    # log-ratio
    zlr <- log2(zz[,"tumor"]*tscl+1) - log2(zz[,"normal"]+1)
    # weights
    # zwts <- log2(1/{1/(zz[,"tumor"]+1) + 1/(zz[,"normal"]+1)})
    zwts <- log2(zz[,"normal"]+1-minBinCount)

    zcna <- CNA(zlr, zchr, zpos, sampleid=sampleid)
    zout <- segment(zcna, weights=zwts)
    zout
}
