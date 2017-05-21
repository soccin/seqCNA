ddir="."
files=dir(ddir,pattern=".R[Dd]ata")
cbs=NULL

computeSegZscore<-function(out) {
  Zscore=rep(NA,nrow(out$output))

  for(ii in seq(nrow(out$output))) {

    idx=out$segRows[ii,1]:out$segRows[ii,2]
    segSd=sd(out$data[idx,3])
    Zscore[ii]=(out$output$seg.mean[ii]/segSd)

  }
  return(Zscore)
}

for(file in files) {
  print(file)
  load(file.path(ddir,file))
  d=out
  d$output[,1]=gsub("___FIXED.*","",d$output[,1])
  d$output$absZscore=abs(round(computeSegZscore(out),2))
  if(is.null(cbs)) {
    cbs=d$output
  } else {
    cbs=rbind(cbs,d$output)
  }
}

cbs$loc.end=as.integer(1e6*cbs$loc.end)
cbs$loc.start=as.integer(1e6*cbs$loc.start)

library(xlsx)

write.xlsx(cbs,
            file=paste0("segTable_",DATE(),".xlsx"),
            row.names=F)
