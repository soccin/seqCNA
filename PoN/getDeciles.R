args=commandArgs(trailing=T)
if(len(args)!=1) {
    cat("\n\tusage: getDeciles.R picardInsStatsDir\n\n")
    quit()
}

SDIR=commandArgs(trailing=F)
SDIR=grep("--file=",SDIR,value=T)
if(len(SDIR)==0) {
    SDIR="."
} else {
    SDIR=dirname(gsub("--file=","",SDIR))
}

require(tidyverse)
source(file.path(SDIR,"read_picard_insstats.R"))
xx=fs::dir_ls(args[1],recurs=T,regex="*.txt") %>% map(read_picard_insStats)
dec=map(xx,get_deciles) %>%
    bind_rows() %>%
    mutate(Decile=seq(11)) %>%
    gather(Sample,Value,-Decile) %>%
    mutate(Sample=basename(Sample) %>% gsub("___INS.*","",.)) %>%
    mutate(Sample=gsub("\\.rg\\.md\\.abra.*","",Sample))
write_csv(dec,"deciles.csv")
