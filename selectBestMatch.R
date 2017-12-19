library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(magrittr)

SNAME=Sys.getenv("SNAME")
SDIR=Sys.getenv("SDIR")
SDIR=ifelse(SDIR=="","seqCNA",SDIR)
source(file.path(SDIR,"include/tools.R"))

args=commandArgs(trailing=T)
if(len(args)!=1) {
    cat("usage: root of seqCNA output directories\n")
    quit()
}

pCon=pipe(paste("find",args[1],"|fgrep .out"))
outFiles=scan(pCon,"")
close(pCon)

zz=list()
ii=1

for(outfile in outFiles) {
    yy=parseSeqCNAOutFile(outfile)
    if(!is.null(yy)) {
        zz[[ii]]=yy
        ii=ii+1
        print(c(ii,outfile))
    }
}

tbl=bind_rows(zz)
tbl %<>% select(which(apply(tbl,2,function(x){len(unique(x))})>1))
tbl %<>% type_convert

bestPairs = tbl %>%
    mutate(TumorId=gsub("__.*","",sampleId)) %>%
    group_by(TumorId) %>% top_n(1,desc(rms.logr)) %>%
    ungroup %>%
    select(sampleId)

bestMatchDirs = tbl %>%
    filter(sampleId %in% bestPairs$sampleId) %>%
    mutate(bestMatchDir=dirname(countFile)) %>%
    select(bestMatchDir) %$%
    as.vector(bestMatchDir)

write(bestMatchDirs,paste0("bestMatches____",gsub("/","___",args[1])))
