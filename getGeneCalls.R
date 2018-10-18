#
# getGeneCalls.R
#

VERSION="5.dev"

SNAME=Sys.getenv("SNAME")
SDIR=Sys.getenv("SDIR")
SDIR=ifelse(SDIR=="","seqCNA",SDIR)
source(file.path(SDIR,"include/parseArgs.R"))
source(file.path(SDIR,"include/tools.R"))
source(file.path(SDIR,"include/globals.R"))


args=list(
    ASSAY=NULL,
    INPUTS=NULL,
    ODIR="."
    )

cat("###############################################################\n")
cat("#\n# Version:",VERSION,"\n")

args=parseArgs(args)

if(is.null(args)) {

    cat(readLines(file.path(SDIR,"docs",paste0(SNAME,".doc"))),sep="\n")
    cat("\n\n")
    quit()

}

#############################################################################
library(data.table)

if(args$ODIR!="." & !dir.exists(args$ODIR)) {
    dir.create(args$ODIR,recursive=T)
}


outDirs=scan(args$INPUTS,"")

# Get Genome from one outfile
out1file=dir(outDirs[1],pattern="_seqSeg.out",full.names=T)
out1=parseSeqCNAOutFile(out1file)

geneAnnoteFile=file.path(
    SDIR,
    glb$geneAnnotationsDir,
    args$ASSAY,
    normalizeGenomeTag(out1$genome),
    "gene_annotations.txt.gz")

geneDb=fread(paste0("gzcat ", geneAnnoteFile))
setkey(geneDb,chrom,start,stop)

getSegFileFromOutDir <- function(odir) {
    segFile=dir(odir,pattern="_seqSeg.seg",full.names=T)
}

segs=rbindlist(lapply(lapply(outDirs,getSegFileFromOutDir),fread))
segs$chrom=as.character(segs$chrom)
setkey(segs,chrom,loc.start,loc.end)

ff=foverlaps(geneDb,segs)

# Get rid of overlaps that do not intersect segments
ff=ff[!is.na(ID)]


#
# DEL priority over DIPLOID over AMP
# ie partial DEL is called but not partial AMP
#

ff=ff[order(ff$seg.mean),]
ff=ff[!duplicated(ff,by=c("ID","gene"))]

require(dtplyr)
require(dplyr)
require(tidyr)
require(tibble)
require(readr)


as_tibble(ff) %>%
    mutate(Zscore=(seg.mean)) %>%
    select(ID,chrom,loc.start,loc.end,gene,transcript,Zscore) %>%
    filter(abs(Zscore)>=.5) %>%
    spread(ID,Zscore,fill="") %>%
    arrange(chrom,loc.start) -> geneEvents

write_csv(geneEvents,file.path(args$ODIR,paste0("SegmentMatrix.csv")))

robustMin <- function(xx) {

    x1=as.numeric((xx[!is.na(xx) & xx!=""]))
    ifelse(len(x1)>0,min(x1),0)

}

iso1=read_tsv("/opt/common/CentOS_6-dev/vcf2maf/v1.6.12/data/isoform_overrides_uniprot")
iso2=read_tsv("/opt/common/CentOS_6-dev/vcf2maf/v1.6.12/data/isoform_overrides_at_mskcc")
canonicalIsoforms=union(iso1[[1]],iso2[[1]])

geneCalls <- geneEvents %>%
    filter(transcript %in% canonicalIsoforms) %>%
    group_by(gene,transcript) %>%
    summarize_at(vars(matches("^s_")),robustMin)

write_csv(geneCalls,file.path(args$ODIR,paste0("GeneMatrix.csv")))

