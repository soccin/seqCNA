#
# checkRlib.R -- check R and the packages that set the output against
# Rlib.lock.tsv
#
# Usage (bin/rscript.sh runs it before every R script):
#
#    Rscript --vanilla bin/checkRlib.R [--list] LOCK LIBDIR
#
# Rlib.lock.tsv (package, version, role) lists only R and the packages
# that compute what ends up in the output; the role column says what
# each computes. The other packages the scripts load (RECORD below)
# can only make a run fail, not change its numbers, and are not
# checked. Stops with exit 1 if R or a locked package has another
# version, if seqDNAcopy does not load from LIBDIR, or if LIBDIR holds
# anything but seqDNAcopy. With --list, also prints one line per locked
# and RECORD package (name, version, where it loads from) for
# Rlib/VERSIONS.txt.
#
# The lock is edited by hand, at a planned version scan only.
#

#
# Loaded by the scripts, recorded in Rlib/VERSIONS.txt, not checked
#
RECORD=c("IRanges", "Rcpp", "RJSONIO", "stringr", "data.table", "dtplyr",
         "dplyr", "tidyr", "tibble", "readr", "openxlsx", "readxl", "rlang", "fs",
         "purrr", "yaml")

fatal <- function(...) {
    cat("\n    FATAL ERROR: ",...,"\n\n",sep="",file=stderr())
    quit(status=1)
}

args=commandArgs(trailingOnly=TRUE)
listVersions=FALSE
if(length(args)>0 && args[1]=="--list") {
    listVersions=TRUE
    args=args[-1]
}
if(length(args)!=2) {
    cat("\n   usage: Rscript --vanilla bin/checkRlib.R [--list] LOCK LIBDIR\n\n",file=stderr())
    quit(status=1)
}
if(!file.exists(args[1])) {
    fatal("no lock file [",args[1],"]")
}
lock=read.delim(args[1],colClasses="character")
libDir=args[2]
if(!dir.exists(libDir)) {
    fatal("no R library [",libDir,"]; run 00.SETUP.sh in the seqCNA checkout")
}
libDir=normalizePath(libDir)

pkgVersion <- function(p) {
    path=suppressWarnings(find.package(p,quiet=TRUE))
    if(length(path)==0) {
        return(c(version=NA,where=NA))
    }
    c(version=unname(read.dcf(file.path(path,"DESCRIPTION"),fields="Version")[1,1]),
      where=dirname(normalizePath(path)))
}

problems=character(0)

rWant=lock$version[lock$package=="R"]
rHave=as.character(getRversion())
if(rHave!=rWant) {
    problems=c(problems,paste0("R ",rHave," (",R.home(),"); locked: ",rWant))
}

extra=setdiff(list.dirs(libDir,full.names=FALSE,recursive=FALSE),"seqDNAcopy")
if(length(extra)>0) {
    problems=c(problems,paste0(libDir," holds more than seqDNAcopy: ",paste(extra,collapse=" ")))
}

pkgs=lock[lock$package!="R",]
for(ii in seq_len(nrow(pkgs))) {
    p=pkgs$package[ii]
    have=pkgVersion(p)
    if(is.na(have["version"])) {
        problems=c(problems,paste0(p,": not installed; locked: ",pkgs$version[ii]))
    } else if(have["version"]!=pkgs$version[ii]) {
        problems=c(problems,paste0(p,": ",have["version"]," in ",have["where"],"; locked: ",pkgs$version[ii]))
    } else if(p=="seqDNAcopy" && have["where"]!=libDir) {
        problems=c(problems,paste0(p,": loads from ",have["where"],", not ",libDir))
    }
}

if(length(problems)>0) {
    cat("\n    FATAL ERROR: R or a locked package does not match ",args[1],":\n",sep="",file=stderr())
    cat(paste0("        ",problems,"\n"),sep="",file=stderr())
    cat("\n    Load R ",rWant," and install the locked versions; seqDNAcopy: run 00.SETUP.sh\n\n",sep="",file=stderr())
    quit(status=1)
}

if(listVersions) {
    cat("R",rHave,R.home(),"locked\n")
    for(p in pkgs$package) {
        have=pkgVersion(p)
        cat(p,have["version"],have["where"],"locked\n")
    }
    for(p in setdiff(RECORD,pkgs$package)) {
        have=pkgVersion(p)
        cat(p,have["version"],have["where"],"\n")
    }
}
