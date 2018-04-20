library(tidyverse)
library(data.table)
library(magrittr)
library(stringr)

canonical=read_tsv("ucsc_mm10_knownCanonical.txt.gz")
canonicalTranscriptsUCSC=canonical %>% distinct(transcript) %>% pull
canonicalTranscripts=read_tsv("ucsc_mm10_knownToEnsembl.txt.gz") %>%
    rename(ID=`#name`) %>%
    filter(ID %in% canonicalTranscriptsUCSC) %>%
    distinct(value) %>%
    pull


GTFFILE="/ifs/depot/annotation/M.musculus/ensembl/v83/Mus_musculus.GRCm38.83.gtf"
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
    rename(chrom=X1,start=X4,stop=X5,gene=gene_name,transcript=transcript_id,exon=exon_number)

write_tsv(iList,"gene_annotations.txt")

