ageSubet = function(ageRange=c(5000,4000),ageMatrix=ageMatrix,ageRef=ages,sites=sites)
{
  agr = as.numeric(colnames(ageMatrix))
  agr.index = which(agr<=ageRange[1]&agr>=ageRange[2])
  ageMatrix.sub = ageMatrix[,agr.index]
  ageMatrix.sub.rowsum=apply(ageMatrix.sub,1,sum)
  ageRef=ageRef[which(ageMatrix.sub.rowsum>0)]
  res = subset(sites,Chronology%in%ageRef)
  return(res)
}