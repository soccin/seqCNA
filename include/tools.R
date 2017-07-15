cc <- function(...) {paste(...,sep='_')}

getGITTag <- function(SDIR) {
    scan(
        pipe(
            paste0(
                "git --git-dir=",
                SDIR,
                "/.git --work-tree=",
                SDIR,
                " describe --tags --always --long --dirty='-UNCOMMITED'")
            ),
        "")
}
