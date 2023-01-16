args=commandArgs(trailing=T)
if(len(args)<1) {
    cat("\n\tusage: findPoN01.R resultsDir\n\n")
    quit()
}

require(tidyverse)
source("seqCNA/include/tools.R")

badPairs=scan("badPairs","")
badNormals=scan("badNormals1","") %>% gsub("/","",.)

xx=fs::dir_ls(args[1],recur=T,regex="\\.out$") %>% map(parseSeqCNAOutFile) %>% map(as_tibble) %>% bind_rows %>% type_convert

qL.MAD=quantile(xx$global.mad,0.025)
qH.MAD=quantile(xx$global.mad,0.975)

# qL.dn=quantile(xx$rms.derivative.noise,0.025)
# qH.dn=quantile(xx$rms.derivative.noise,0.975)

qL.nb=quantile(xx$numBins,0.025)
qH.nb=quantile(xx$numBins,0.975)

xf=xx %>%
    filter(global.mad>qL.MAD & global.mad<qH.MAD) %>%
    filter(numBins>qL.nb & numBins<qH.nb)

ds=xf %>% select(sampleId,numBins,rms.logr,numSegments,wRMSD.seg.mean.AUTO) %>% separate(sampleId,c("Nt","Nc"),sep="__") %>% mutate(Pt=gsub("_N[0-9]+_d.*","",Nt)) %>% mutate(Pc=gsub("_N[0-9]+_d.*","",Nc)) %>% filter(Pt!=Pc) %>% select(-Pt,-Pc) %>% gather(NType,SID,Nt,Nc) %>% arrange(wRMSD.seg.mean.AUTO)

dx=ds %>%
    mutate(R1=rank(wRMSD.seg.mean.AUTO)) %>%
    group_by(SID) %>%
    summarize(
        med.wRMSD=median(wRMSD.seg.mean.AUTO),
        med.nSeg=median(numSegments),
        N=n()) %>%
    arrange(med.wRMSD) %>%
    arrange(desc(N))

pon=dx %>% mutate(BadNorm=SID %in% badNormals) %>% filter(med.wRMSD<.07 & !BadNorm)

pon1=pon %>% pull(SID)

write(pon1,"pon1")

test.samples=xx %>% select(sampleId,numBins,rms.logr,numSegments,wRMSD.seg.mean.AUTO) %>% separate(sampleId,c("Nt","Nc"),sep="__") %>% filter(Nt %in% pon1 & Nc %in% pon1) %>% unite(sampleId,Nt,Nc,sep="__") %>% pull(sampleId)

write(test.samples,"test.samples")
