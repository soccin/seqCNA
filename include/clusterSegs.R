##########################################
#
# Cluster segments
# taken from facets facets-clustersegs.R (commit: f501272)
#
# DMP uses threshold=0.08
# 2^(-2*0.04) ~ 0.95

clustersegs <- function(out,threshold=0.08) {

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
        ocnlevels <- ocnlevels[-j]
        ocnclust[ocnclust>j] <- ocnclust[ocnclust>j] - 1
        segid[segid > j] <- segid[segid > j] - 1
        wj=weights[segid==j]/sum(weights[segid==j])
        ocnlevels[j] <- sum(wj*cnlr[segid==j])
    }


    nProbesPerSeg=table(segid)
    diploidClusterNum=as.numeric(which.min(abs(ocnlevels/sqrt(nProbesPerSeg))))

    out$cluster=list(
        diploidClusterNum=diploidClusterNum,
        originalLevels=ocnlevels)

    ocnlevels=ocnlevels - ocnlevels[diploidClusterNum]

    segs$cluster=ocnclust
    segs$clust.mean=ocnlevels[segs$cluster]
    segs=segs[order(segs$segNo),]

    deltaDat=rep(segs$clust.mean-segs$seg.mean,segs$num.mark)

    out$output$seg.mean=segs$clust.mean
    out$output$cluster=segs$cluster
    out$data[,3]=out$data[,3]+deltaDat

    out
}
