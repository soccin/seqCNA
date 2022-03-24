STATES=3

#colors5State=RColorBrewer::brewer.pal(5,"RdBu")

gg=MetBrewer::scale_color_met_d("OKeeffe1")
colors5State=gg$palette(11)[c(1,4,6,8,11)]

if(STATES==5) {
    cnvStates=(c(2,1,0,-1,-2))
    cnvLevels=c(0,-1,-2,1,2)
    cnvColors=(colors5State)
} else {
    cnvStates=(c(1,0,-1))
    cnvLevels=c(0,-1,1)
    cnvColors=colors5State[c(1,3,5)]
}

require(tidyverse)
require(readxl)
require(fs)
require(patchwork)
library(gridExtra)

if(STATES==3) {
    relevelCNV<-function(x){factor(sign(x),levels=cnvLevels)}
} else {
    relevelCNV<-function(x){factor(x,levels=cnvLevels)}
}

getEventTbl<-function(dx,ti) {
    filter(dx,TYPE==ti) %>%
    select(Gene,Sample,CNV) %>%
    group_by(Gene) %>%
    summarize(
        pAmp=paste0(round(100*mean(CNV>0),1),"%"),nAmp=sum(CNV>0),
        pDel=paste0(round(100*mean(CNV<0),1),"%"),nDel=sum(CNV<0)
        ) %>%
    arrange(desc(Gene))
}

cArgs=commandArgs(trailing=T)

dd=read_xlsx(cArgs[1])

if(STATES==3) {
    dd=dd %>% mutate(CNV=sign(seg.mean))
} else {
    rlang::abort("Not implemented 5 state")
}

nSamps=dd %>% distinct(ID) %>% nrow

genes=dd %>%
    group_by(gene) %>%
    summarize(nEvents=n(),PCT=nEvents/nSamps) %>%
    arrange(desc(PCT)) %>%
    filter(PCT>.1) %>%
    pull(gene)

dd=dd %>% select(ID,CNV,Gene=gene) %>% filter(Gene %in% genes)

xx=dd %>% spread(Gene,CNV,fill=0)

geneOrder=genes

subTypeOrder=rep(0,nSamps)

xs=xx

for(ig in seq(geneOrder)) {
    print(ig)
    subTypeOrder=5*subTypeOrder+as.numeric(relevelCNV(xs[[geneOrder[ig]]]))-1
}

samp.order=xs$ID[order(subTypeOrder,decreasing=T)]

halt("INCLUDE")

dx=dd %>% mutate(fCNV=factor(CNV,levels=cnvStates),ID=factor(ID,levels=samp.order),Gene=factor(Gene,levels=rev(geneOrder))) %>% filter(Gene %in% genes)
pg1=dx %>% ggplot(aes(ID,Gene)) + theme_light(base_size=16) + geom_tile(aes(fill=fCNV),color="grey50") + scale_fill_manual(values=cnvColors,drop=F) +  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())

pdf(file=cc("argosCNVOncoPrintV1.pdf"),width=11,height=8.5)
print(pg1)
dev.off()
