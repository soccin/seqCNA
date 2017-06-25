################################################################
#
# This code will parse command line args in the form of
#    KEY=VAL
# and sets
#    args[[KEY]]=VAL
#
# If an arg is required set it to NULL. And if any are missing
# the function will print a message and return NULL.
#
# All args are converted to UPPERCASE
#

suppressPackageStartupMessages(require(stringr))

parseArgs<-function(args,convertCaseArgs=T) {

    cArgs=commandArgs(trailing=T)

    if(length(cArgs)>0) {

        parseArgs=str_match(cArgs,"(.*)=(.*)")

        if(convertCaseArgs) {
            parseArgs[,2]=toupper(parseArgs[,2])
        }

        dummy=apply(parseArgs,1,function(x){args[[str_trim(x[2])]]<<-str_trim(x[3])})

    }

    if(any(sapply(args,is.null))) {
        missing=which(sapply(args,is.null))
        cat("\nmissing required arg(s)\n\n    ")
        for(ii in missing){
            cat(names(args)[[ii]],"")
        }
        cat("\n\n")
        return(NULL)
    }

    return(args)

}

