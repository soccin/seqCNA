fixSampleNames<-function(x) {

    x=gsub("indelRealigned_recal_","",x)
    x=gsub("___MD","",x)
    x=gsub("^s_","",x)
    x=gsub(".bam$","",x)
    return(x)

}
