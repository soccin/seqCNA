suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
})

compute_deciles <- function(bam) {

    sid=get_sm_tag(bam)

    insFile=picard_CollectInsertSizeMetrics(bam)
    stats=read_picard_insStats(insFile)
    dec=tibble(extract_deciles(stats))
    colnames(dec)=sid

    dec

}


picard_CollectInsertSizeMetrics<-function(bam) {

    oFile=tempfile()
    eFile=tempfile()
    res=system2(
        "picardV2",
        c(
            "CollectInsertSizeMetrics",
            paste0("I=",bam),
            "H=/dev/null",
            paste0("O=",oFile)
        ),
        stderr=ifelse(interactive(),"",eFile)
    )

    if(res>0) {

        cat("\n\tFATAL ERROR: compute_deciles::picard_CollectInsertSizeMetrics\n")
        cat("\tRES =",res,"\n")
        cat("\tBAM =",bam,"\n")
        cat("\n")
        fs::file_copy(eFile,cc("log","compute_deciles","picard","CISM",DATE(),".stderr"))
        rlang::abort("FATAL ERROR")

    }

    oFile

}

get_sm_tag<-function(bam) {
    system2(file.path(SDIR,"getSampleTag.sh"),bam,stdout=T)
}

extract_deciles <- function(insStats) {

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
    pcols=cols(
      insert_size = col_double(),
    All_Reads.fr_count = col_double()
    )

    skip=grep("^## HISTOGRAM",readLines(pfile,20))
    read_tsv(pfile,skip=skip,col_types=pcols,progress=F) %>%
        filter(All_Reads.fr_count>0)
}


if(!exists("SDIR")) {
    SDIR=Sys.getenv("SDIR")
    if(SDIR=="") {
        SDIR="seqCNA"
    }
}