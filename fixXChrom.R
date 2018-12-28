suppressPackageStartupMessages(library(tidyverse))
args=commandArgs(trailingOnly=T)
seg=read_tsv(args[1])
maxChrom=seg %>% distinct(chrom) %>% max
cat("\n\nMax Chromosome =",maxChrom,"==> X\n\n")
seg.new=seg %>% mutate(chrom=ifelse(chrom==maxChrom,"X",as.character(chrom)))
write_tsv(seg.new,args[1])
