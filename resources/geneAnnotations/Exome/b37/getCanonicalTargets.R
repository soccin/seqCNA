library(tidyverse)
library(data.table)
library(magrittr)
library(stringr)


ISOFORM0="/opt/common/CentOS_6-dev/vcf2maf/v1.6.12/data/isoform_overrides_at_mskcc"
ISOFORM1="/opt/common/CentOS_6-dev/vcf2maf/v1.6.12/data/isoform_overrides_uniprot"

isoform0=read_tsv(ISOFORM0,
    col_names=c("TID","gene_name","refseq_id","ccds_id"),skip=1)
isoform1=read_tsv(ISOFORM1,
    col_names=c("TID","gene_name","refseq_id","ccds_id"),skip=1)

isoforms=bind_rows(isoform0,isoform1) %>% distinct(gene_name,.keep_all=T)
canonicalTranscripts = isoforms %>% distinct(TID) %>% pull

# canonical=read_tsv("ucsc_mm10_knownCanonical.txt.gz")
# canonicalTranscriptsUCSC=canonical %>% distinct(transcript) %>% pull
# canonicalTranscripts=read_tsv("ucsc_mm10_knownToEnsembl.txt.gz") %>%
#     rename(ID=`#name`) %>%
#     filter(ID %in% canonicalTranscriptsUCSC) %>%
#     distinct(value) %>%
#     pull


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
    rename(chrom=X1,start=X4,stop=X5,gene=gene_name,transcript=transcript_id,exon=exon_number)

ozfile=gzfile("gene_annotations.txt.gz","w",compression=9)
write_tsv(iList,ozfile)
close(ozfile)

