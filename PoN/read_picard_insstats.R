get_deciles <- function(insStats) {

    cs=cumsum(insStats$All_Reads.fr_count)

    dec.vec = rep(NA, 11)
    dec.vec[1] = min(insStats$insert_size)
    dec.vec[11] = max(insStats$insert_size)

    for(i in 1:9){
        nth = (i*(sum(insStats$All_Reads.fr_count) + 1))/10
        idx = which(cs < round(nth))
        dec.vec[i+1] = insStats$insert_size[rev(idx)[1] + 1]
    }

    dec.vec

}

read_picard_insStats <- function(pfile) {
    skip=grep("^## HISTOGRAM",readLines(pfile,20))
    read_tsv(pfile,skip=skip) %>%
        filter(All_Reads.fr_count>0)
}

compute_decile_ssdev <- function(dm){
    dt=dm %>% spread(Sample,Value)
    samps=colnames(dt)[-1]
    nSamps=len(samps)
    SS=list()
    for(i in seq(nSamps)) {
        for(j in seq(nSamps)) {
            if(i<j) {
                SS[[len(SS)+1]]=tibble(Si=samps[i],Sj=samps[j],ssDev=sum((dt[,i+1]-dt[,j+1])^2))
            }
        }
    }
    bind_rows(SS)
}