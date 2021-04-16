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

    #browser()
    #
    # This was the old way of finding the shift to try
    # to set the diploid level. Failed if there was a tiny
    # segment very close to zero even there was a cluster
    # of many more segments(points) far enough away.
    # Ie if one segment of a few points just happend to have a
    # mean of zero then there would be no shift since
    #     smallButNonZeroMean/sqrt(VeryLargeNumberOfPoints) > 0
    #
    # nProbesPerSeg=table(segid)
    # diploidClusterNum=as.numeric(which.min(abs(ocnlevels/sqrt(nProbesPerSeg))))

    # out$cluster=list(
    #     diploidClusterNum=diploidClusterNum,
    #     originalLevels=ocnlevels)

    # ocnlevels=ocnlevels - ocnlevels[diploidClusterNum]
    # }

    #
    # New method minimize
    # $$
    #    \sum_i \{(s_i-\Delta)n_i\}^2
    # $$
    # for $\delta$. This is:
    # $$
    #    \Delta=\frac{\sum_i s_i n^2_i}{\sum_i  n^2_i}
    # $$
    #
    # n_i == nProbesPerSeg
    # s_i == ocnlevels
    #
    # But one more complexity; if we have a genome with a large number of events
    # more then the number of points that diploid we do not to shift that event to
    # so we only look at segments that are within one MAD of the scatter of
    # probes to their means

    #
    # MAD of (probes-seg.mean)
    #
    data.mad=mad(out$data[,3]-rep(segs$seg.mean,segs$num.mark))
    ocn.ii=which(abs(ocnlevels)<data.mad)

    nProbesPerSeg=table(segid)
    Delta=sum(ocnlevels[ocn.ii]*(nProbesPerSeg[ocn.ii])^2)/sum(nProbesPerSeg[ocn.ii]^2)

    ocnlevels=ocnlevels - Delta

    segs$cluster=ocnclust
    segs$clust.mean=ocnlevels[segs$cluster]
    segs=segs[order(segs$segNo),]

    deltaDat=rep(segs$clust.mean-segs$seg.mean,segs$num.mark)

    out$output$seg.mean=segs$clust.mean
    out$output$cluster=segs$cluster
    out$data[,3]=out$data[,3]+deltaDat

    out
}
