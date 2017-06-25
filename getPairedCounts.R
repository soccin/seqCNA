#
# getPairedCounts.R
#

VERSION="5.dev"

SDIR=Sys.getenv("SDIR")
SDIR=ifelse(SDIR=="",".",SDIR)
source(file.path(SDIR,"include/parseArgs.R"))
source(file.path(SDIR,"include/misc.R"))

args=list(
    BINSIZE="auto",
    TUMOR=NULL,
    NORMAL=NULL,
    MAPLOC="false",
    GENOME="hg19",
    ODIR="."
    )

args=parseArgs(args)

cat("###############################################################\n")
cat("#\n# Version:",VERSION,"\n")

if(is.null(args)) {
    cat("\n\nusage: getPairCounts TUMOR=/path/tumor.bam NORMAL=/path/normal.bam BINSIZE=[auto] MAPLOC=[false] GENOME=[hg19]\n\n")
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

tBase=fixSampleNames(basename(tumorBam))
nBase=fixSampleNames(basename(normalBam))

if(tBase==nBase){
    stop("tumor==normal")
}

sampleId=cc(tBase,"_",nBase)
cat("# sampleId =",sampleId,"\n")

bb=bams2counts(normalBam,tumorBam,X=T,gbuild=args$GENOME)
