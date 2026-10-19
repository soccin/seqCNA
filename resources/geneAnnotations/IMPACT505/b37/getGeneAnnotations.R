# https://cmo.mskcc.org/cmo/wp-content/uploads/2020/07/IMPACT505_Gene_list_detailed.xlsx

require(tidyverse)
require(readxl)

genes505=read_xlsx("IMPACT505_Gene_list_detailed.xlsx",sheet=2) %>% rename_all(~gsub(" ","_",.))

geneSyn=read_csv("badNames.csv")

genes505=left_join(genes505,select(geneSyn,Approved_Symbol,GeneCards),by="Approved_Symbol") %>%
    mutate(GeneSym=ifelse(is.na(GeneCards),Approved_Symbol,GeneCards)) %>%
    arrange(GeneCards,GeneSym)

ISOFORM0="/opt/common/CentOS_6-dev/vcf2maf/v1.6.12/data/isoform_overrides_at_mskcc"
ISOFORM1="/opt/common/CentOS_6-dev/vcf2maf/v1.6.12/data/isoform_overrides_uniprot"

isoform0=read_tsv(ISOFORM0,
    col_names=c("TID","gene_name","refseq_id","ccds_id"),skip=1)
isoform1=read_tsv(ISOFORM1,
    col_names=c("TID","gene_name","refseq_id","ccds_id"),skip=1)

isoforms=bind_rows(isoform0,isoform1) %>% distinct(gene_name,.keep_all=T)

# left_join(genes505,isoforms,by=c(Approved_Symbol="gene_name")) %>%
#     filter(is.na(TID)) %>%
#     select(1:4) %>%
#     write_csv("badNames.csv")

genes505=left_join(genes505,isoforms,by=c(GeneSym="gene_name"))

canonicalTranscripts = genes505 %>% distinct %>% pull(TID)

GTFFILE="/ifs/depot/annotation/H.sapiens/ensembl/v75/Homo_sapiens.GRCh37.75.gtf"
gtf=read_tsv(GTFFILE,col_names=F,comment="#",col_types=list(X1 = col_character()))

iList=gtf %>%
    mutate(transcript_id=str_match(X9,'transcript_id "([^;]*)";')[,2]) %>%
    mutate(gene_name=str_match(X9,'gene_name "([^;]*)";')[,2]) %>%
    mutate(exon_number=str_match(X9,'exon_number "(\\d+)";')[,2]) %>%
    select(-X9) %>%
    filter(X3=="exon") %>%
    filter(transcript_id %in% canonicalTranscripts) %>%
    mutate(X1=gsub("^chr","",X1)) %>%
    select(X1,X4,X5,gene_name,transcript_id,exon_number) %>%
    rename(chrom=X1,start=X4,stop=X5,gene=gene_name,transcript=transcript_id,exon=exon_number) %>%
    left_join(isoforms,by=c(transcript="TID")) %>%
    mutate(transcript=gsub("\\.\\d+$","",refseq_id),exon=paste0("exon",exon)) %>%
    select(1:6) %>%
    mutate(chrom=as.numeric(ifelse(chrom=="X",23,chrom))) %>%
    arrange(chrom,start)

write_tsv(iList,"gene_annotations.txt.gz")
