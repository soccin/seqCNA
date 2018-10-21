library(tidyverse)
impactIlist="/home/socci/Work/CMO/Mouse_IMPACT/BICTargetFiles/M-IMPACT_v1_mm10/M-IMPACT_v1_mm10_targets.ilist"

header=readLines(impactIlist,1000)
header=grep("^@",header,value=T)
if(len(header)==1000) {
    stop("Need bigger cache")
}

ann <- read_tsv(impactIlist,col_names=F,comment="@") %>%
    filter(grepl(":ENSMUST",X5)) %>%
    separate(X5,c("gene","transcript"),sep=":") %>%
    arrange(transcript,X2) %>%
    group_by(gene,transcript) %>%
    mutate(exon=paste0("exon",row_number())) %>%
    rename(chrom=X1,start=X2,stop=X3) %>%
    select(chrom,start,stop,gene,transcript,exon) %>%
    ungroup %>%
    mutate(chrom=factor(gsub("^chr","",chrom),levels=c(1:99,"X","Y","M","MT"))) %>%
    arrange(chrom,start)

write_tsv(ann,"M-IMPACT_v1_mm10.txt")
