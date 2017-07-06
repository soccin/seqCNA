#
# getPairedCounts.R
#

VERSION="5.dev"

SNAME=Sys.getenv("SNAME")
SDIR=Sys.getenv("SDIR")
SDIR=ifelse(SDIR=="",".",SDIR)
source(file.path(SDIR,"include/parseArgs.R"))
source(file.path(SDIR,"include/misc.R"))

args=list(
    TUMOR=NULL,
    NORMAL=NULL,
    GENOME="hg19",
    ODIR=".",
    GCNORM=TRUE,
    SAMPLEID=""
    )


cat("###############################################################\n")
cat("#\n# Version:",VERSION,"\n")

args=parseArgs(args)

if(is.null(args)) {

    cat(readLines(file.path(SDIR,"docs",paste0(SNAME,".doc"))),sep="\n")
    cat("\n\n")
    quit()

}

if(args$ODIR!="." & !dir.exists(args$ODIR)) {
    dir.create(args$ODIR,recursive=T)
}

################################################################
# Load rest of libraries after args checkout
#

suppressPackageStartupMessages(library(seqDNAcopy))
cat("# Version(seqDNAcopy):",sessionInfo()$otherPkgs$seqDNAcopy$Version,"\n")
keys=sort(names(args))
for(key in keys) {
    cat("#",key,"=",args[[key]],"\n")
}

################################################################

normalBam=args$NORMAL
tumorBam=args$TUMOR

if(args$SAMPLEID=="") {

    tBase=fixSampleNames(basename(tumorBam))
    nBase=fixSampleNames(basename(normalBam))

    if(tBase==nBase){
        stop("tumor==normal")
    }

    sampleId=paste(tBase,"_",nBase,sep="_")

} else {

    sampleId=args$SAMPLEID

}

cat("# sampleId =",sampleId,"\n")

bb=bams2counts(normalBam,tumorBam,X=T,gbuild=args$GENOME,GCcorrect=args$GCNORM)

d=list()
d$sampleId=sampleId
d$bb=bb
d$args=args
dat=list()
dat[[make.names(sampleId)]]=d

save(dat,file=file.path(args$ODIR,paste0(sampleId,".rda")),compress=T)


