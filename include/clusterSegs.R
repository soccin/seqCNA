##########################################
#
# Cluster segments
# taken from facets facets-clustersegs.R (commit: f501272)
#
# DMP uses threshold=0.1
# 2^(-2*0.04) ~ 0.95

clustersegs <- function(out,threshold=0.04) {

    segs=out$output
    nsegs=nrow(segs)
    segs$segNo=seq(nsegs)

    cnlr=out$data[is.finite(out$data[,3]),3]
    weights=out$weights[is.finite(out$data[,3])]

    segid <- rep(rank(segs$seg.mean, ties.method="random"), segs$num.mark)

    cnlr <- cnlr[order(segid)]
    weights <- weights[order(segid)]
    segid <- sort(segid)

    segs=segs[order(segs$seg.mean),]

    ocnclust<-1:nsegs
    ocnlevels<-segs$seg.mean

    while ((length(ocnlevels) > 1) && (min(diff(ocnlevels)) < threshold)) {
        j <- which.min(diff(ocnlevels))
        print(j)
        print(which(ocnlevels==1.0710))
        ocnlevels <- ocnlevels[-j]
        ocnclust[ocnclust>j] <- ocnclust[ocnclust>j] - 1
        segid[segid > j] <- segid[segid > j] - 1
        wj=weights[segid==j]/sum(weights[segid==j])
        ocnlevels[j] <- mean(wj*cnlr[segid==j])

        print(ocnlevels)
    }

    #out$output$cluster=ocnclust[rank(out$output$seg.mean,ties.method="random")]
    #out$output$clust.mean=round(ocnlevels[out$output$cluster],4)
