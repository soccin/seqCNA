#
# installRlib.R -- install the packages of Rlib.lock.tsv into a library
#
# Usage (00.SETUP.sh runs it):
#
#    Rscript --vanilla bin/installRlib.R LOCK LIBDIR CACHEDIR [OLDCACHE]
#
# Each package is installed from its source tarball, at the locked
# version, in lock order; a tarball whose md5 differs from the lock
# stops the install. Tarballs are kept in CACHEDIR; one already in
# OLDCACHE (the src/ of the previous Rlib/) is copied, not downloaded.
# The R version must match the lock, and LIBDIR must be the only library
# besides R's own (00.SETUP.sh sets R_LIBS, R_LIBS_USER and R_LIBS_SITE
# to it), so nothing is taken from another library. The submodule row
# (seqDNAcopy) is left to 00.SETUP.sh.
#

scriptFile=normalizePath(sub("^--file=","",grep("^--file=",commandArgs(),value=TRUE)))
source(file.path(dirname(scriptFile),"rlibTools.R"))

args=commandArgs(trailingOnly=TRUE)
if(!(length(args) %in% 3:4)) {
    cat("\n   usage: Rscript --vanilla bin/installRlib.R LOCK LIBDIR CACHEDIR [OLDCACHE]\n\n",file=stderr())
    quit(status=1)
}
lock=readLock(args[1])
libDir=normalizePath(args[2])
cacheDir=args[3]
oldCache=if(length(args)==4 && dir.exists(args[4])) args[4] else NULL
dir.create(cacheDir,showWarnings=FALSE,recursive=TRUE)

rWant=lock$version[lock$source=="R"]
if(as.character(getRversion())!=rWant) {
    rlibFatal("R ",as.character(getRversion())," on PATH; the lock needs R ",rWant)
}

if(!identical(normalizePath(.libPaths()),normalizePath(c(libDir,.Library)))) {
    rlibFatal("library paths are not [",libDir,"] and R's own: ",paste(.libPaths(),collapse=" "))
}

pkgs=lock[!(lock$source=="R" | startsWith(lock$source,"submodule:")),]

for(ii in seq_len(nrow(pkgs))) {
    row=pkgs[ii,]
    tarball=fetchTarball(row,cacheDir,oldCache)
    cat(sprintf("[%d/%d] %s %s (%s)\n",ii,nrow(pkgs),row$package,row$version,row$source))
    install.packages(tarball,lib=libDir,repos=NULL,type="source",quiet=TRUE)
    desc=file.path(libDir,row$package,"DESCRIPTION")
    installed=if(file.exists(desc)) read.dcf(desc,fields="Version")[1,1] else "none"
    if(installed!=row$version) {
        rlibFatal("install of ",row$package," ",row$version," failed (installed: ",installed,")")
    }
}
