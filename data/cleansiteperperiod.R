cleaning=read.csv("data/siteperpref.csv")

remK <- function(i)gsub(",","",i)
sinK <- apply(cleaning[,2:ncol(cleaning)],2,remK)
cleaning[,2:ncol(cleaning)]=apply(sinK,2,as.numeric)
write.csv(file="data/siteperprefperiod_cleaned.csv",cleaning,row.names=F)
cleaned=read.csv("data/siteperperprefperiod_cleaned.csv")
