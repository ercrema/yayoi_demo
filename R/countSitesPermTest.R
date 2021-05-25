countSitePermTest = function(ageRange=c(8000,2000),ageMatrix=ageMatrix,ageRef=ages,sites=sites,nsim=100,ncores=2,backsight=100,regions)
{
  require(doSNOW)
  require(foreach)
  require(Matrix.utils)
  
  each.region = sites$Region
  K = length(regions)
  tt = length(ageRange[1]:ageRange[2])
  combinedAgeMatrix = ageMatrix[match(sites$Chronology,ageRef),which(as.numeric(colnames(ageMatrix))<=ageRange[1]&as.numeric(colnames(ageMatrix))>=ageRange[2])]
  siteIDs = unique(sites$SiteID)
  rownames(combinedAgeMatrix)=sites$SiteID
  observed = matrix(NA,ncol=K,nrow=length(ageRange[1]:ageRange[2]))
  
  print('Computing Observed Site Counts per Region...')
  for (i in 1:K)
  {
    sitesID.region = unique(sites$SiteID[which(each.region==regions[i])])
    combinedAgeMatrix.region = combinedAgeMatrix[which(rownames(combinedAgeMatrix)%in%sitesID.region),]
    agMat.region=aggregate.Matrix(combinedAgeMatrix.region, row.names(combinedAgeMatrix.region))
    agMat.region=agMat.region>0
    observed[,i]=apply(agMat.region,2,sum)
  }
  

  print('Permuting Sites...')
  pb <- txtProgressBar(min=0, max=nsim, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  cl <- snow::makeCluster(ncores)
  registerDoSNOW(cl)
  on.exit(snow::stopCluster(cl))	
  
  
  sim.counts = foreach (i=1:nsim,.options.snow = opts,.packages=c('Matrix.utils'),.combine='c') %dopar% {
    sites.region.tmp = unique(data.frame(SiteID=sites$SiteID,Region=sites$Region))
    sites.region.tmp$Region = sample(sites.region.tmp$Region)
    sites.random.region=merge(sites,sites.region.tmp,by.x='SiteID',by.y='SiteID',all.x=TRUE)
    sim.counts.regions = matrix(NA,ncol=K,nrow=length(ageRange[1]:ageRange[2]))
    
    for (j in 1:K)
    {
      sitesID.region = unique(sites.random.region$SiteID[which(sites.random.region$Region.y==regions[j])])
      combinedAgeMatrix.region = combinedAgeMatrix[which(rownames(combinedAgeMatrix)%in%sitesID.region),]
      agMat.region=aggregate.Matrix(combinedAgeMatrix.region, row.names(combinedAgeMatrix.region))
      agMat.region=agMat.region>0
      sim.counts.regions[,j]=apply(agMat.region,2,sum)
    }
 return(sim.counts.regions)
  }

  sim.counts.mat=array(sim.counts,dim=c(tt,K,nsim))
  
  #Convert to annual growth rate: expression((t1/t0)^(1/d)-1)
  timeSequence=ageRange[1]:ageRange[2]
  foo = function(x,backsight,timeSequence)
  {
    obs=rep(NA,length(timeSequence))
    for (i in 1:c(length(obs)-backsight))
    {
      d=backsight 	
      t0 = x[i]
      t1 = x[i+backsight]
      obs[i+backsight] = (t1/t0)^(1/d)-1
      if (t1==0|t0==0){obs[i+backsight]=NA}
    }
    return(obs)
  }
  
  observed.growth = apply(observed,2,foo,timeSequence=timeSequence,backsight=backsight)

  results = vector('list',length=K)
  results.growth = vector('list',length=K)
  
  
  for (i in 1:K)
  {
    sim.growth = apply(sim.counts.mat[,i,],2,foo,timeSequence=timeSequence,backsight=backsight)
    sim.cnts = sim.counts.mat[,i,]
    
    lo.growth = apply(sim.growth,1,quantile,0.025,na.rm=T)
    hi.growth = apply(sim.growth,1,quantile,0.975,na.rm=T)
    lo.cnts = apply(sim.cnts,1,quantile,0.025,na.rm=T)
    hi.cnts = apply(sim.cnts,1,quantile,0.975,na.rm=T)
    
    mean.growth=apply(sim.growth,1,mean)
    sd.growth=apply(sim.growth,1,sd)
    z.growth =(sim.growth-mean.growth)/sd.growth
    mean.cnts=apply(sim.cnts,1,mean)
    sd.cnts=apply(sim.cnts,1,sd)
    z.cnts =(sim.cnts-mean.cnts)/sd.cnts
    
    z.lo.cnts = apply(z.cnts,1,quantile,0.025,na.rm=T)
    z.hi.cnts = apply(z.cnts,1,quantile,0.975,na.rm=T)
    observed.z.cnts=(observed[,i]-mean.cnts)/sd.cnts
    z.lo.growth = apply(z.growth,1,quantile,0.025,na.rm=T)
    z.hi.growth = apply(z.growth,1,quantile,0.975,na.rm=T)
    observed.z.growth=(observed.growth[,i]-mean.growth)/sd.growth
    
    

    observed.sumstat.cnts = sum((observed.z.cnts-z.hi.cnts)[which(observed.z.cnts>z.hi.cnts)])+sum((z.lo.cnts-observed.z.cnts)[which(observed.z.cnts<z.lo.cnts)])
    sim.sumstat.cnts=apply(z.cnts,2,function(x,z.lo,z.hi){return(sum((x-z.hi)[which(x>z.hi)])+sum((z.lo-x)[which(x<z.lo)]))},z.lo=z.lo.cnts,z.hi=z.hi.cnts)
    observed.sumstat.growth = sum((observed.z.growth-z.hi.growth)[which(observed.z.growth>z.hi.growth)])+sum((z.lo.growth-observed.z.growth)[which(observed.z.growth<z.lo.growth)])
    sim.sumstat.growth=apply(z.growth,2,function(x,z.lo,z.hi){return(sum((x-z.hi)[which(x>z.hi)])+sum((z.lo-x)[which(x<z.lo)]))},z.lo=z.lo.growth,z.hi=z.hi.growth)
    
    
    
    obs.growthrate = data.frame(BP=timeSequence,obs=observed.growth[,i],lo=lo.growth,hi=hi.growth)
    obs.cnts = data.frame(BP=timeSequence,obs=observed[,i],lo=lo.cnts,hi=hi.cnts)
    
    r.cnts = sum(observed.sumstat.cnts<=sim.sumstat.cnts)
    r.growth = sum(observed.sumstat.growth<=sim.sumstat.growth)
    
    pval.growth = (r.growth+1)/(nsim+1)
    pval.cnts = (r.cnts+1)/(nsim+1)
    
    results[[i]]=list(obs.growthrate=obs.growthrate,obs.cnts=obs.cnts,pval.growth=pval.growth,pval.cnts=pval.cnts,sim.cnts=sim.cnts)
  }
  
  return(results)
}