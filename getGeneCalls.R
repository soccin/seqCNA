#
# getGeneCalls.R
#

convertGeneSymbolsMouseToHuman <- function(mgg) {

    source(file.path(SDIR,"convert_mouse_to_human.R"))

    mgg=tibble(Gene=unique(mgg))
    genesV2=convert_mouse_to_human(mgg) %>%
        rename(HGNC.symbol=Gene) %>%
        bind_cols(mgg) %>%
        rename(MGI.symbol=Gene)

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

    cat(readLines(file.path(SDIR,"docs","usage",paste0(SNAME,".doc"))),sep="\n")
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

genome=normalizeGenomeTag(out1$genome)

#
# Gene calls are made on the autosomes only. A call on X or Y needs a
# sex-matched normal and the known sex of the sample, which we do not
# have. The segments name X by number (23 human, 20 mouse) and the
# annotations name it X, 23 or not at all, so keep the autosomes by
# number in both rather than rely on the names not matching.
#

autosomes=switch(genome,
    b37=as.character(1:22),
    mm10=as.character(1:19),
    {
        cat("\n\n    FATAL ERROR: getGeneCalls: unknown genome [",genome,"]\n\n",sep="")
        quit(status=1)
    })

geneAnnoteFile=file.path(
    SDIR,
    glb$geneAnnotationsDir,
    args$ASSAY,
    genome,
    "gene_annotations.txt.gz")

geneDb=fread(cmd=paste0("zcat ", geneAnnoteFile))
geneDb$chrom=as.character(geneDb$chrom)
geneDb=geneDb[chrom %in% autosomes]

setkey(geneDb,chrom,start,stop)

getSegFileFromOutDir <- function(odir) {
    segFile=dir(odir,pattern="_seqSeg.seg",full.names=T)
}

segs=rbindlist(lapply(lapply(outDirs,getSegFileFromOutDir),fread))
segs$chrom=as.character(segs$chrom)
segs=segs[chrom %in% autosomes]
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
    cat("\n    No significant copy number changes\n\n",file=readmeFile)
    cat("\n    No significant copy number changes\n\n")
    quit()
}


write_csv(geneEvents,file.path(args$ODIR,paste0("SegmentMatrix.csv")),na="")

robustMin <- function(xx) {

    x1=as.numeric((xx[!is.na(xx) & xx!=""]))
    ifelse(len(x1)>0,min(x1),0)

}

#
# Canonical isoforms for Exome: the vcf2maf v1.6.12 override lists,
# vendored in resources/vcf2maf/v1.6.12/ (they were read from
# /opt/common on the LSF clusters).
#
if(args$ASSAY %in% c("Exome")) {
    isoDir=file.path(SDIR,"resources/vcf2maf/v1.6.12")
    iso1=read_tsv(file.path(isoDir,"isoform_overrides_uniprot"),show_col_types=FALSE)
    iso2=read_tsv(file.path(isoDir,"isoform_overrides_at_mskcc"),show_col_types=FALSE)
    canonicalIsoforms=union(iso1[[1]],iso2[[1]])
} else {
    canonicalIsoforms=sort(unique(geneEvents$transcript))
}

geneCalls <- geneEvents %>%
    filter(transcript %in% canonicalIsoforms) %>%
    group_by(gene,transcript) %>%
    summarize_at(vars(matches("^s_")),robustMin)

if(args$ASSAY=="M-IMPACT_v2") {
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


if(args$ASSAY=="M-IMPACT_v2") {
    geneTable=geneTable %>%
        left_join(mm.genes,by=c(gene="MGI.symbol")) %>%
        select(ID,gene,HGNC.gene=HGNC.symbol,everything())
}

require(openxlsx)
pID=basename(gsub("/$","",args$ODIR))
write.xlsx(geneTable,file.path(args$ODIR,paste0(pID,"___","GeneTable.xlsx")))
