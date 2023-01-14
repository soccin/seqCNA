cc <- function(...) {paste(...,sep='_')}

args=commandArgs(trailing=T)
if(length(args)!=3) {
    cat("\n\tusage: getBestDecileNormals.R DECILE_DB BAM WORK_DIR\n\n")
    quit()
}

library(dplyr)
library(tidyr)
library(tibble)

SNAME=Sys.getenv("SNAME")
SDIR=Sys.getenv("SDIR")
SDIR=ifelse(SDIR=="","seqCNA",SDIR)

source(file.path(SDIR,"include/compute_deciles.R"))

decileDb=args[1]
bam=args[2]
WDIR=args[3]

deciles=read_csv(decileDb)

dec=compute_deciles(bam)

dM=deciles %>% filter(Sample!=names(dec)) %>% spread(Decile,Value) %>% data.frame %>% column_to_rownames("Sample") %>% t

bestNormals=names(sort(apply(dM,2,function(x){sum((x-dec[[1]])^2)})) %>% head(20))
write(bestNormals,file.path(WDIR,cc("bestNormals","_",get_sm_tag(bam))))
