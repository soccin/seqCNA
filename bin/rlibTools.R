#
# rlibTools.R -- helpers shared by mkRlibLock.R, installRlib.R and
# checkRlib.R. Base R only: these run before Rlib/ exists.
#

rlibFatal <- function(...) {
    cat("\n    FATAL ERROR: ",...,"\n\n",sep="",file=stderr())
    quit(status=1)
}

#
# Rlib.lock.tsv: package, version, source, md5. source is R (the R
# version itself), CRAN, Bioc/<release>, github:<owner>/<repo>@<commit>
# or submodule:<dir>.
#
readLock <- function(lockFile) {
    if(!file.exists(lockFile)) {
        rlibFatal("no lock file [",lockFile,"]")
    }
    read.delim(lockFile,colClasses="character")
}

#
# Where the source tarball of one lock row can be downloaded. A CRAN
# package moves to Archive/ when a newer version is released.
#
tarballUrls <- function(row) {
    tb=paste0(row$package,"_",row$version,".tar.gz")
    if(row$source=="CRAN") {
        cran="https://cloud.r-project.org/src/contrib"
        c(file.path(cran,tb),file.path(cran,"Archive",row$package,tb))
    } else if(startsWith(row$source,"Bioc/")) {
        bioc=file.path("https://bioconductor.org/packages",sub("^Bioc/","",row$source),"bioc/src/contrib")
        c(file.path(bioc,tb),file.path(bioc,"Archive",row$package,tb))
    } else if(startsWith(row$source,"github:")) {
        repo=sub("@.*$","",sub("^github:","",row$source))
        commit=sub("^.*@","",row$source)
        paste0("https://github.com/",repo,"/archive/",commit,".tar.gz")
    } else {
        rlibFatal("no download for source [",row$source,"] of ",row$package)
    }
}

#
# Return the path of the row's tarball in cacheDir, downloading it if
# it is not there. With row$md5 set, a cached or downloaded file must
# match it. oldCache, if given, is searched before downloading.
#
fetchTarball <- function(row,cacheDir,oldCache=NULL) {
    dest=file.path(cacheDir,paste0(row$package,"_",row$version,".tar.gz"))
    md5ok <- function(f) {
        file.exists(f) && (is.na(row$md5) || unname(tools::md5sum(f))==row$md5)
    }
    if(md5ok(dest)) {
        return(dest)
    }
    if(!is.null(oldCache) && md5ok(file.path(oldCache,basename(dest)))) {
        file.copy(file.path(oldCache,basename(dest)),dest,overwrite=TRUE)
        return(dest)
    }
    for(url in tarballUrls(row)) {
        rc=tryCatch(download.file(url,dest,mode="wb",quiet=TRUE),
                    error=function(e) 1,warning=function(w) 1)
        if(rc==0 && md5ok(dest)) {
            return(dest)
        }
        if(rc==0) {
            rlibFatal(row$package," ",row$version,": md5 of ",url," is not ",row$md5)
        }
    }
    unlink(dest)
    rlibFatal("cannot download ",row$package," ",row$version," from: ",paste(tarballUrls(row),collapse=" "))
}
