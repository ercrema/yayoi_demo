library(ggplot2)
library(nimbleCarbon)
library(rcarbon)
library(rnaturalearth)
library(maptools)
library(here)

# Setup Window of Analysis and Site Locations ----

## Generate Window
sf_japan <- ne_states(country = "japan") |> subset(name_vi != "Okinawa")
window <- as(sf_japan, "SpatialPolygons") |>  unionSpatialPolygons(IDs = rep(1, nrow(sf_japan))) |> spTransform(CRSobj = CRS("+proj=utm +zone=54 ellps=WGS84"))

## Generate Spatial Units
source(here('R','gridder.R'))
grd  <- gridder(window,cell_diameter=120000,method='hex')
# plot(grd);length(grd)

## Compute Distance Matrix
grd.centers  <-  gCentroid(grd,byid = TRUE)
distMatrix  <- coordinates(grd.centers) |> dist() |> as.matrix()
distMatrix  <- distMatrix/1000 

## Assign Sites to Grids (NOT RELEVANT HERE)
# locations  <- spsample(n=1000,type='random',x=window)
# mat  <- t(gIntersects(locations,grd,TRUE))
# index  <- as.numeric(unlist(apply(mat,1,function(x){which(x==TRUE)})))
# plot(locations,pch=1,col='lightgrey')
# points(locations[which(index==as.numeric(names(sort(table(index),decreasing=T)[1])))],pch=20)
# locations <- as(locations,'SpatialPointsDataFrame')
# locations@data = data.frame(SiteID=1:length(locations),GridID=index)

# Simulate Data ----
## Generate Samples Sizes for Each Region
avgSample = 1000
nSamples <- rpois(length(grd),lambda=avgSample/length(grd))
id.grid  <- rep(1:length(grd),nSamples)

## Covariance Function
cov_GPL2 <- nimbleFunction(
  run = function(dists = double(2), rhosq = double(0), etasq = double(0), sigmasq = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    deltaij <- matrix(nrow = n, ncol = n,init = TRUE)
    diag(deltaij) <- 1
    for(i in 1:n)
      for(j in 1:n)
        result[i, j] <- etasq*exp(-rhosq*dists[i,j]^2)+sigmasq*deltaij[i,j]
    return(result)
  })
Ccov_GPL2 <- compileNimble(cov_GPL2)

## Spatially Structured Growth Model
spatial_demo  <- nimbleCode({
	for (i in 1:nGrids){
		local_r[i]  <- r + s[i]
	}
	for (i in 1:nDates){
		theta[i] ~ dExponentialGrowth(a=start,b=end,r=local_r[id.grid[i]])
	}
	mu_s[1:nGrids] <- 0;
	cov_s[1:nGrids, 1:nGrids] <- cov_GPL2(dmat[1:nGrids, 1:nGrids], rhosq, etasq, 0.0000001)
	s[1:nGrids] ~ dmnorm(mu_s[1:nGrids], cov = cov_s[1:nGrids, 1:nGrids])
})

### Parameters
etasq <-0.00001
rhosq <- 0.0001
r <- 0.004
curve(etasq * exp(-rhosq * x^2), from = 0, to = 200)


### Actual Simulation
set.seed(1234)
simModel <- nimbleModel(code = spatial_demo,constants = list(start=3000,end=1800,nDates = sum(nSamples),dmat = distMatrix,etasq=etasq,rhosq=rhosq,r=r,nGrids=length(grd),id.grid=id.grid))
simModel$simulate('s')
simModel$simulate('local_r')
simModel$simulate('theta')

### Aggregate and Generate Final DF
d  <- data.frame(ID=1:sum(nSamples),GridID=id.grid,cra=round(uncalibrate(round(simModel$theta))$ccCRA),cra_error=20)
d$medDates  <- medCal(calibrate(d$cra,d$cra_error))

### Sanity Checks
grd@data$r = r+simModel$s
grd@data$id = 1:length(grd)
grd.sf  <-  as(grd,'sf')
grd.sf  <- cbind(grd.sf,st_coordinates(st_centroid(grd.sf)))
ggplot(grd.sf,aes(fill=r))  + geom_sf() + scale_fill_gradient2()

maxRegion=which.max(grd@data$r)
minRegion=which.min(grd@data$r)

maxRegionSPD <- calibrate(d$cra[which(d$GridID==grd@data$id[maxRegion])],d$cra_error[which(d$GridID==grd@data$id[maxRegion])]) |> spd(timeRange=c(3000,1800))
minRegionSPD <- calibrate(d$cra[which(d$GridID==grd@data$id[minRegion])],d$cra_error[which(d$GridID==grd@data$id[minRegion])]) |> spd(timeRange=c(3000,1800))
par(mfrow=c(2,1))
plot(maxRegionSPD)
plot(minRegionSPD)






# Analyse Simulated Data ----

### Define Constants
data(intcal20)
constants <- list(N=sum(nSamples),Ngrids=length(grd),dmat=distMatrix,calBP=intcal20$CalBP,C14BP=intcal20$C14Age,C14err=intcal20$C14Age.sigma,a=3000,b=1800,GridID=id.grid,a=3000,b=1800)

### Define Data
data <- list(cra=d$cra,cra_error=d$cra_error)

### Define Init
inits  <- list(s=rep(0,constants$Ngrids),etasq=0.0001,rhosq=0.001,r=0.0004,theta=d$medDates)
inits$cov_s <- Ccov_GPL2(constants$dists, inits$rhosq, inits$etasq, 0.0000001)
inits$s <-  t(chol(inits$cov_s)) %*% rnorm(constants$Ngrids)
inits$s <- inits$s[ , 1]  # so can give nimble a vector rather than one-column matrix
inits$local_r <- inits$r + inits$s
inits$theta[which(inits$theta>constants$a)]  <- constants$a
inits$theta[which(inits$theta<constants$b)]  <- constants$b

### Define Core Model
spatial_model  <- nimbleCode({
	for (i in 1:Ngrids){
		local_r[i]  <- r + s[i]/10;
	}
	for (i in 1:N){
		theta[i] ~ dExponentialGrowth(a=a,b=b,r=local_r[GridID[i]])
		c14age[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
		sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
		sigmaDate[i] <- (cra_error[i]^2+sigmaCurve[i]^2)^(1/2);
		cra[i] ~ dnorm(mean=c14age[i],sd=sigmaDate[i]);
	}
	mu_s[1:Ngrids] <- 0;
	cov_s[1:Ngrids, 1:Ngrids] <- cov_GPL2(dmat[1:Ngrids, 1:Ngrids], rhosq, etasq, 0.0000001)
	s[1:Ngrids] ~ dmnorm(mu_s[1:Ngrids], cov = cov_s[1:Ngrids, 1:Ngrids])
	r ~ dnorm(0,sd=0.001)
	etasq ~ dexp(10000)
	rhosq ~ dexp(10000)
})

## MCMC (Single Core Test) ----

model <- nimbleModel(spatial_model,constants = constants,data=data,inits=inits)
cModel <- compileNimble(model)
conf <- configureMCMC(model)
conf$addMonitors('s')
conf$removeSamplers('s[1:85]')
conf$addSampler(c('r','s[1:85]'), type='RW_block') 
MCMC <- buildMCMC(conf)
cMCMC <- compileNimble(MCMC)
testRes <- runMCMC(cMCMC, nchain=1,niter = 5000, thin=1,nburnin = 1000,samplesAsCodaMCMC = T,progressBar = TRUE) 

## MCMC (Multi-Core Parallel) ----






