require(tidyverse)
require(tidygenomics)

ILISTFILE="/juno/projects/BIC/targets/designs/M-IMPACT_v2_mm10/M-IMPACT_v2_mm10_targets.ilist"

geneTargets=read_tsv(ILISTFILE,comment="@",col_names=F) %>%
    filter(!grepl("^rs\\d+",X5)) %>%
    mutate(X1=gsub("chr","",X1)) %>%
    mutate(GID=row_number())
genes=geneTargets %>% distinct(X5) %>% arrange(X5) %>% pull

MM_GTF="/juno/depot/annotation/M.musculus/ensembl/v83/Mus_musculus.GRCm38.83.gtf"

gtf0=read_tsv(MM_GTF,comment="#",col_names=F)

gtf=gtf0 %>%
    filter(X3=="exon") %>%
    select(1,4,5,9) %>%
    mutate(RID=row_number()) %>%
    separate_rows(X9,sep="; ") %>%
    filter(grepl("^gene_name|^transcript_id|^exon_number",X9)) %>%
    separate(X9,c("key","val"),sep=" ") %>%
    mutate(val=gsub('"',"",val)) %>%
    spread(key,val) %>%
    rename(X2=X4,X3=X5)


#    filter(gene_name %in% genes) %>%

oo=genome_intersect(geneTargets,gtf)

o1=oo %>% left_join(geneTargets %>% select(GID,GS=X2,GE=X3)) %>% left_join(gtf %>% select(RID,ES=X2,EE=X3)) %>% rowwise %>% mutate(Overlap=min(EE,GE)-max(ES,GS)) %>% ungroup %>% group_by(GID) %>% slice_max(Overlap)


# > setdiff(genes,g1)
# [1] "Exosc6" "Kmt5a"

canonicalT=o1 %>% group_by(X5)  %>% count(transcript_id  ) %>% arrange(desc(n)) %>% distinct(X5,.keep_all=T) %>% pull(transcript_id)
canonicalAnnote=o1 %>% filter(transcript_id %in% canonicalT)

annote=gtf %>%
    filter(transcript_id %in% canonicalT) %>%
    mutate(exon_number=paste0("exon",exon_number)) %>%
    select(chrom=X1,start=X2,stop=X3,gene=gene_name,transcript=transcript_id,exon=exon_number) %>%
    filter(chrom!="X") %>% type_convert %>% arrange(chrom,start)

write_tsv(annote,"gene_annotations.txt.gz")
