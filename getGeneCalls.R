#
# getGeneCalls.R
#

VERSION="5.dev"

SNAME=Sys.getenv("SNAME")
SDIR=Sys.getenv("SDIR")
SDIR=ifelse(SDIR=="",".",SDIR)
source(file.path(SDIR,"include/parseArgs.R"))

FIX THIS

GET GENOME FROM OUTPUT FILES NORMALIZE HG19 ==> B37

args=list(
    ASSAY=NULL,


require(data.table)
require(xlsx)
require(stringr)

GTF_RDA="gencode.vM13.annotation.sorted.gtf.rda"

if(!file.exists(GTF_RDA)) {
    MM10_GTF="/ifs/depot/annotation/M.musculus/gencode/vM13/gencode.vM13.annotation.sorted.gtf.gz"
    gtfCols=c("chrom","db","feature","start","end","score","strand","frame","attribute")

    gtf=fread(paste("zcat",MM10_GTF),col.names=gtfCols)
    gtf[,chrom:=gsub("chr","",chrom)]
    gtf[,chrom:=as.numeric(factor(chrom,c(1:19,"X","Y","M")))]

    gtfFlds=c("gene_name","transcript_id","exon_number","gene_type","transcript_type")
    for(gfld in gtfFlds) {
        print(gfld)
        gtf[,(gfld):=gsub('("|;)','',str_match(attribute,paste0(gfld," (.*?) "))[,2])]
    }

    gtf_attribute=gtf$attribute

    gtf[,attribute:=NULL]
    gtf[,length:=end-start]

    save(gtf,gtf_attribute,file=GTF_RDA)
} else {
    load(GTF_RDA)
}

gene.gtf=gtf
gene.gtf=gene.gtf[
                    gene_type=="protein_coding"
                    & transcript_type=="protein_coding"
                    & feature %in% "transcript",]
gene.gtf=gene.gtf[!duplicated(gene.gtf$gene_name),]

setkey(gene.gtf,chrom,start,end)

seg=read.xlsx("proj_07593_SegTable.xlsx",sheetIndex=1)
seg=data.table(seg)
setkey(seg,chrom,loc.start,loc.end)

ff=foverlaps(gene.gtf,seg)

# Get rid of overlaps that do not intersect segments
ff=ff[!is.na(ID)]


#
# DEL priority over DIPLOID over AMP
# ie partial DEL is called but not partial AMP
#

ff=ff[order(ff$seg.mean),]
ff=ff[!duplicated(ff,by=c("ID","gene_name"))]

require(dplyr)
require(tidyr)
require(tibble)

as_tibble(ff) %>%
    mutate(Zscore=sign(seg.mean)*absZscore) %>%
    select(ID,chrom,loc.start,loc.end,gene_name,transcript_id,Zscore) %>%
    filter(abs(Zscore)>=1) %>%
    spread(ID,Zscore,fill="") %>%
    filter(!grepl("Rik\\d*$",gene_name)) %>%
    arrange(chrom,loc.start) -> geneEvents

require(xlsx)
write.xlsx2(as.data.frame(geneEvents),"proj_07593__GeneEvents.xlsx",row.names=F)

