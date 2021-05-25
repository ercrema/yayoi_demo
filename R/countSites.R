countSites = function(ageRange=c(10000,2000),ageMatrix=ageMatrix,ageRef=ages,sites=sites)
{
  require(Matrix.utils)
  combinedAgeMatrix = ageMatrix[match(sites$Chronology,ageRef),which(as.numeric(colnames(ageMatrix))<=ageRange[1]&as.numeric(colnames(ageMatrix))>=ageRange[2])]
  siteIDs = unique(sites$SiteID)
  rownames(combinedAgeMatrix)=sites$SiteID
  agMat=aggregate.Matrix(combinedAgeMatrix, row.names(combinedAgeMatrix))
  agMat=agMat>0
  res=data.frame(BP=ageRange[1]:ageRange[2],Cnts=apply(agMat,2,sum))
  return(res)
}