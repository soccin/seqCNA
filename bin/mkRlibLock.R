#
# mkRlibLock.R -- write Rlib.lock.tsv from an installed R library
#
# Usage, from the repo directory, at a planned version scan only:
#
#    Rscript --vanilla bin/mkRlibLock.R LIBDIR CACHEDIR > Rlib.lock.tsv
#
# LIBDIR is the library holding the versions to pin (with R's own
# library behind it); CACHEDIR receives the source tarballs, whose md5
# goes into the lock. The lock lists R itself, then every package the
# pipeline loads with its whole Depends/Imports/LinkingTo closure,
# except R's base packages, in install order; seqDNAcopy (the
# submodule) is last. 00.SETUP.sh installs exactly these into Rlib/ and
# bin/checkRlib.R checks them before every R script.
#
# After writing a new lock: rerun 00.SETUP.sh, rerun a reference
# project, compare with the last one, and record the change in
# CHANGELOG.md and METHODS.md.
#

#
# Packages the pipeline scripts load (library, require, ::) and
# seqDNAcopy imports. Keep in step with the scripts.
#
TOP=c("DNAcopy", "pctGCdata", "Rsamtools", "IRanges", "Rcpp", "RJSONIO",
      "stringr", "data.table", "dtplyr", "dplyr", "tidyr", "tibble", "readr",
      "openxlsx", "readxl", "rlang", "fs")

#
# Packages not installed from CRAN or Bioconductor: owner/repo@commit
#
GITHUB=c(pctGCdata="veseshan/pctGCdata@d2d4fafd10595750977f9bf5ebfeaed6e78da781")

BIOC="3.22"

scriptFile=normalizePath(sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)))
BINDIR=dirname(scriptFile)
SDIR=dirname(BINDIR)
source(file.path(BINDIR,"rlibTools.R"))

args=commandArgs(trailingOnly=TRUE)
if(length(args)!=2) {
    cat("\n   usage: Rscript --vanilla bin/mkRlibLock.R LIBDIR CACHEDIR > Rlib.lock.tsv\n\n",file=stderr())
    quit(status=1)
}
libDir=args[1]
cacheDir=args[2]
dir.create(cacheDir,showWarnings=FALSE,recursive=TRUE)

ip=installed.packages(lib.loc=c(libDir,.Library),fields="biocViews")
ip=ip[!duplicated(ip[,"Package"]),]
base=rownames(installed.packages(lib.loc=.Library,priority="base"))

deps=tools::package_dependencies(TOP,db=ip,which=c("Depends","Imports","LinkingTo"),recursive=TRUE)
pkgs=setdiff(unique(c(TOP,unlist(deps))),base)

missing=setdiff(pkgs,rownames(ip))
if(length(missing)>0) {
    rlibFatal("not installed: ",paste(missing,collapse=" "))
}

#
# Install order: a package after everything it depends on
#
direct=tools::package_dependencies(pkgs,db=ip,which=c("Depends","Imports","LinkingTo"))
order=character(0)
while(length(order)<length(pkgs)) {
    left=setdiff(pkgs,order)
    ready=left[vapply(left,function(p) all(intersect(direct[[p]],pkgs) %in% order),logical(1))]
    if(length(ready)==0) {
        rlibFatal("dependency cycle among: ",paste(left,collapse=" "))
    }
    order=c(order,sort(ready))
}

source=vapply(order,function(p) {
    if(p %in% names(GITHUB)) {
        paste0("github:",GITHUB[[p]])
    } else if(!is.na(ip[p,"biocViews"])) {
        paste0("Bioc/",BIOC)
    } else {
        "CRAN"
    }
},character(1))

lock=data.frame(package=order,version=unname(ip[order,"Version"]),source=unname(source),md5=NA_character_)

for(ii in seq_len(nrow(lock))) {
    tarball=fetchTarball(lock[ii,],cacheDir)
    lock$md5[ii]=unname(tools::md5sum(tarball))
}

seqDNAcopyVersion=read.dcf(file.path(SDIR,"seqDNAcopy","DESCRIPTION"),fields="Version")[1,1]

lock=rbind(
    data.frame(package="R",version=as.character(getRversion()),source="R",md5="-"),
    lock,
    data.frame(package="seqDNAcopy",version=seqDNAcopyVersion,source="submodule:seqDNAcopy",md5="-"))

write.table(lock,stdout(),sep="\t",quote=FALSE,row.names=FALSE)
