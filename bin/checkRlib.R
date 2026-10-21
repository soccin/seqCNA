#
# checkRlib.R -- check R and Rlib/ against Rlib.lock.tsv
#
# Usage (bin/rscript.sh runs it before every R script):
#
#    Rscript --vanilla bin/checkRlib.R [--list] LOCK LIBDIR
#
# Stops with exit 1 if R is not the locked version, or if any locked
# package is missing from LIBDIR, has another version there, or would
# be loaded from another library. With --list, also prints one line per
# package (name, version, source, where it loads from) for
# Rlib/VERSIONS.txt.
#

scriptFile=normalizePath(sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)))
source(file.path(dirname(scriptFile),"rlibTools.R"))

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
lock=readLock(args[1])
libDir=args[2]

problems=character(0)

rWant=lock$version[lock$source=="R"]
rHave=as.character(getRversion())
if(rHave!=rWant) {
    problems=c(problems,paste0("R ",rHave," (",R.home(),"); locked: ",rWant))
}

if(!dir.exists(libDir)) {
    rlibFatal("no R library [",libDir,"]; run 00.SETUP.sh in the seqCNA checkout")
}
libDir=normalizePath(libDir)

pkgs=lock[lock$source!="R",]
where=character(nrow(pkgs))
for(ii in seq_len(nrow(pkgs))) {
    p=pkgs$package[ii]
    path=suppressWarnings(find.package(p,quiet=TRUE))
    if(length(path)==0) {
        problems=c(problems,paste0(p,": not installed; locked: ",pkgs$version[ii]))
        next
    }
    where[ii]=dirname(normalizePath(path))
    have=read.dcf(file.path(path,"DESCRIPTION"),fields="Version")[1,1]
    if(where[ii]!=libDir) {
        problems=c(problems,paste0(p,": loads from ",where[ii],", not ",libDir))
    } else if(have!=pkgs$version[ii]) {
        problems=c(problems,paste0(p,": ",have,"; locked: ",pkgs$version[ii]))
    }
}

if(length(problems)>0) {
    cat("\n    FATAL ERROR: R or its packages do not match ",args[1],":\n",sep="",file=stderr())
    cat(paste0("        ",problems,"\n"),sep="",file=stderr())
    cat("\n    R must be ",rWant," (module load R/4.5.1_nds); rebuild the library with 00.SETUP.sh\n\n",sep="",file=stderr())
    quit(status=1)
}

if(listVersions) {
    cat("R",rHave,R.home(),"\n")
    for(ii in seq_len(nrow(pkgs))) {
        cat(pkgs$package[ii],pkgs$version[ii],pkgs$source[ii],where[ii],"\n")
    }
}
