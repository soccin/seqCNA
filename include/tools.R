cc <- function(...) {paste(...,sep='_')}

getGITTag <- function(SDIR) {

    pCon=pipe(
        paste0(
            "git --git-dir=",
            SDIR,
            "/.git --work-tree=",
            SDIR,
            " describe --tags --always --long --dirty='-UNCOMMITED'"
            )
        )

    gitStr=scan(pCon,"")
    close(pCon)

    gitStr

}

# Parses the
#    key = value
# format of the current seqCNA out file
# returns a list

parseSeqCNAOutFile <- function(outfile) {

    xx=readLines(outfile)
    if(any(grepl("X.seg.mean.avg",xx))) {

        yy=sapply(
            lapply(xx[-1],function(x){strsplit(x," = ")[[1]]}),
            function(x){ll=list();ll[[x[1]]]=gsub(" $","",x[2]);ll})

        return(yy)

    } else {

        return(NULL)

    }

}

# Normalize Genome
# Need a better fix but for now since
# we are using hg19/b37 interchangably
# normalize to b37
#

normalizeGenomeTag <- function(gg) {
    ifelse(gg=="hg19","b37",gg)
}


