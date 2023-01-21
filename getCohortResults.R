args=commandArgs(trailing=T)
if(len(args)<1) {
    cat("\n   usage: getCohortGeneTable.R ProjectOutputFolder\n\n")
    quit()
}

require(tidyverse)
require(readxl)

gfiles=fs::dir_ls(args[1],recur=T,regex="Gene.*xlsx$")
geneTbl=map(gfiles,read_xlsx) %>% bind_rows %>% type_convert %>% filter(chrom<=22) %>% separate(ID,c("Tumor","Normal"),sep="__") %>% select(-Normal) %>% arrange(Tumor,chrom,loc.start)
projNo=strsplit(gfiles[1],"/")[[1]][1]

dir.create(file.path(projNo,"COHORT"),recursive=T,showWarnings=F)

openxlsx::write.xlsx(geneTbl,file.path(projNo,"COHORT",cc(projNo,"_","GeneTable_Vx_1_1.xlsx")))

segTbl=fs::dir_ls(args[1],recur=T,regex="IGV.seg$") %>% map(read_tsv) %>% bind_rows %>% mutate(ID=gsub("__s_.*","",ID))

write_tsv(segTbl,file.path(projNo,"COHORT",cc(projNo,"_","Vx_1_1__IGV.seg")))
