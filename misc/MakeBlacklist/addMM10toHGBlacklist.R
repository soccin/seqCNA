args=commandArgs(trailing=T)

# load human blacklists
load(args[1])
humanBLs=ls()[grep("bl$",ls())]

mm10bl=list()

#
# Dummy blacklist with no entries
# conform to VS spec

MAXPOS=2^28
mm=matrix(c(MAXPOS,MAXPOS),nrow=1)
colnames(mm)=c("start","end")

#
# Mouse: auto=1-19, X=20, Y=21, MT=22
#
for(chromi in c(1:19,20,21,22)) {
    print(chromi)
    mm10bl[[chromi]]=mm
}

save(list=c(humanBLs,"mm10bl"),file="hgBlacklist.rda")