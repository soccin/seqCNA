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
#
# Keep the columns that differ between pairs, and always the three used
# below: with a single pair every column is constant.
#
tbl %<>% select(sampleId, countFile, rms.logr, where(~ len(unique(.x)) > 1))
tbl %<>% type_convert

bestPairs = tbl %>%
    mutate(TumorId=gsub("__s_.*","",sampleId)) %>%
    mutate(NormalId=gsub(".*__s_","s_",sampleId)) %>%
    filter(TumorId!=NormalId) %>%
    group_by(TumorId) %>%
    top_n(1,desc(rms.logr)) %>%
    ungroup %>%
    select(sampleId)

bestMatchDirs = tbl %>%
    filter(sampleId %in% bestPairs$sampleId) %>%
    mutate(bestMatchDir=dirname(countFile)) %>%
    select(bestMatchDir) %$%
    as.vector(bestMatchDir)

write(bestMatchDirs,paste0("bestMatches____",gsub("/","___",args[1])))
