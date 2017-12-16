#####
#
# getSDIR():
#   returns the current directory a script is actually in
#   to allow relative imports. No way to boot strap this so
#   need to have it explicitly copied to all scripts that need
#   this functionality, cannot source it.
#
#   Plus it is hideous
#

getSDIR <- function(){
    args=commandArgs(trailing=F)
    TAG="--file="
    path_idx=grep(TAG,args)
    SDIR=dirname(substr(args[path_idx],nchar(TAG)+1,nchar(args[path_idx])))
    if(length(SDIR)==0) {
        return(getwd())
    } else {
        return(SDIR)
    }
}
