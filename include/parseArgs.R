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

        for(ii in seq(nrow(parseArgs))) {

            key=str_trim(parseArgs[ii,2])
            if(convertCaseArgs) {
                key=toupper(key)
            }
            value=str_trim(parseArgs[ii,3])

            if(!is.null(args[[key]])) {

                # If the arg is already set to a default value
                # preserve that type
                args[[key]]=do.call(paste0("as.",class(args[[key]])),list(value))

            } else {

                args[[key]]=value

            }

        }

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

