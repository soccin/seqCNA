#
# selectBestMatch.R
#
# usage: selectBestMatch OUTDIR
#
# For each tumor keep the normal whose pair has the lowest rms.logr.
# The tumor and normal of each pair come from OUTDIR/pipeline_info/pairs.tsv
# (written by PIPE.sh), not from parsing the pair id; rms.logr comes from
# each pair's _seqSeg.out. Writes, in the current directory:
#
#   bestMatches____<OUTDIR>  directory of each best pair, one per line
#   pairing.tsv              TUMOR, NORMAL, PAIR, rms.logr; one row per
#                            tumor. Internal: it names the normals.
#

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
})

SNAME=Sys.getenv("SNAME")
SDIR=Sys.getenv("SDIR")
SDIR=ifelse(SDIR=="","seqCNA",SDIR)
source(file.path(SDIR,"include/tools.R"))

args=commandArgs(trailing=T)
if(len(args)!=1) {
    cat("usage: selectBestMatch OUTDIR (root of seqCNA output directories)\n")
    quit()
}
outRoot=args[1]

pairsFile=file.path(outRoot,"pipeline_info","pairs.tsv")
if(!file.exists(pairsFile)) {
    cat("\n\tFATAL ERROR no pairs file [",pairsFile,"] (written by PIPE.sh)\n\n")
    quit(status=1)
}
pairs=read_tsv(pairsFile,col_types=cols(.default="c"))

pCon=pipe(paste("find",outRoot,"|fgrep .out"))
outFiles=scan(pCon,"")
close(pCon)

tbl=outFiles |>
    lapply(parseSeqCNAOutFile) |>
    bind_rows() |>
    select(sampleId,countFile,rms.logr) |>
    mutate(rms.logr=as.numeric(rms.logr)) |>
    left_join(pairs |> select(PAIR,TUMOR,NORMAL),by=c(sampleId="PAIR"))

unknown=tbl |> filter(is.na(TUMOR)) |> pull(sampleId)
if(len(unknown)>0) {
    cat("\n\tFATAL ERROR pairs not in",pairsFile,":",unknown,"\n\n")
    quit(status=1)
}

#
# One best normal per tumor; a tie in rms.logr goes to the first pair id
#
best=tbl |>
    arrange(sampleId) |>
    group_by(TUMOR) |>
    slice_min(rms.logr,n=1,with_ties=FALSE) |>
    ungroup() |>
    arrange(TUMOR)

write(dirname(best$countFile),paste0("bestMatches____",gsub("/","___",outRoot)))

pairing=best |> select(TUMOR,NORMAL,PAIR=sampleId,rms.logr)
write_tsv(pairing,"pairing.tsv")
print(as.data.frame(pairing))
