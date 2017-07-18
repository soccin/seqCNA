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
