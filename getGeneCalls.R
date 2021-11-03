#
# getGeneCalls.R
#

convertGeneSymbolsMouseToHuman <- function(mgg) {

    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

    mggu=unique(sort(mgg))

    genesV2 = getLDS(
        attributes = c("mgi_symbol"),
        filters = "mgi_symbol",
        values = mggu,
        mart = mouse,
        attributesL = c("hgnc_symbol"),
        martL = human,
        uniqueRows=T)

    genesV2=data.table(genesV2)
    setkey(genesV2,MGI.symbol)

    genesV2

}

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

geneDb=fread(cmd=paste0("zcat ", geneAnnoteFile))
geneDb$chrom=as.character(geneDb$chrom)

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

suppressPackageStartupMessages({
    require("biomaRt")
    require(dtplyr)
    require(dplyr)
    require(tidyr)
    require(tibble)
    require(readr)
})

chromOrder=c(seq(1:99),"X","Y","M","MT")

QCUT=0.05
LOGR_CUT=1

geneEvents <- as_tibble(ff) %>%
    filter(FDR<QCUT & abs(seg.mean)>LOGR_CUT) %>%
    mutate(seg.mean=round(seg.mean,3)) %>%
    spread(ID,seg.mean,fill=NA) %>%
    mutate(chrom=factor(chrom,levels=chromOrder)) %>%
    arrange(chrom,loc.start)

if(nrow(geneEvents)==0) {
    readmeFile=file.path(args$ODIR,"README.txt")
    cat("\n    No siginificant CopyNumber Changes\n\n",file=readmeFile)
    cat("\n    No siginificant CopyNumber Changes\n\n")
    quit()
}


write_csv(geneEvents,file.path(args$ODIR,paste0("SegmentMatrix.csv")),na="")

robustMin <- function(xx) {

    x1=as.numeric((xx[!is.na(xx) & xx!=""]))
    ifelse(len(x1)>0,min(x1),0)

}

if(args$ASSAY %in% c("Exome")) {
    iso1=read_tsv("/opt/common/CentOS_6-dev/vcf2maf/v1.6.12/data/isoform_overrides_uniprot")
    iso2=read_tsv("/opt/common/CentOS_6-dev/vcf2maf/v1.6.12/data/isoform_overrides_at_mskcc")
    canonicalIsoforms=union(iso1[[1]],iso2[[1]])
} else {
    canonicalIsoforms=sort(unique(geneEvents$transcript))
}

geneCalls <- geneEvents %>%
    filter(transcript %in% canonicalIsoforms) %>%
    group_by(gene,transcript) %>%
    summarize_at(vars(matches("^s_")),robustMin)

if(args$ASSAY=="M-IMPACT_v1") {
    mm.genes=convertGeneSymbolsMouseToHuman(geneCalls$gene) %>%
        tibble %>%
        distinct(MGI.symbol,.keep_all=T)

    geneCalls=geneCalls %>%
        left_join(mm.genes,by=c(gene="MGI.symbol")) %>%
        select(gene,HGNC.gene=HGNC.symbol,everything())
}

write_csv(geneCalls,file.path(args$ODIR,paste0("GeneMatrix.csv")))

geneTable <- as_tibble(ff) %>%
    filter(transcript %in% canonicalIsoforms) %>%
    filter(FDR<QCUT & abs(seg.mean)>LOGR_CUT) %>%
    arrange(seg.mean) %>%
    distinct(ID,gene,transcript,.keep_all=T) %>%
    select(ID,gene,transcript,chrom,loc.start,loc.end,seg.mean,FDR) %>%
    mutate(chrom=factor(chrom,levels=chromOrder)) %>%
    mutate(seg.mean=round(seg.mean,3)) %>%
    arrange(ID,chrom,loc.start)


if(args$ASSAY=="M-IMPACT_v1") {
    geneTable=geneTable %>%
        left_join(mm.genes,by=c(gene="MGI.symbol")) %>%
        select(ID,gene,HGNC.gene=HGNC.symbol,everything())
}

require(openxlsx)
pID=basename(gsub("/$","",args$ODIR))
write.xlsx(geneTable,file.path(args$ODIR,paste0(pID,"___","GeneTable.xlsx")))
