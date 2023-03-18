library(ggplot2)
library(nimbleCarbon)
library(rcarbon)
library(rnaturalearth)
library(maptools)
# Setup Window of Analysis and Site Locations ----
sf_japan <- ne_states(country = "japan") |> subset(name_vi != "Okinawa")
window <- as(sf_japan, "SpatialPolygons") |>  unionSpatialPolygons(IDs = rep(1, nrow(sf_japan))) |> spTransform(CRSobj = CRS("+proj=utm +zone=54 ellps=WGS84"))
# Sample some random locations
N  <-  1000
locations <- spsample(window,n=N,type='random')
# Compute Distance Matrix
distMat <- coordinates(locations) |> dist() |> as.matrix()
distMat <- distMat / 1000 #Convert in Km

# Simulate Observed Data----
## Covariance Function ----
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

## Core Model Simulation ----
spatial_demo  <- nimbleCode({
	for (i in 1:N){
		local_r[i]  <- r + s[i]/10
		theta[i] ~ dExponentialGrowth(a=start,b=end,r=local_r[i])
	}
	mu_s[1:N] <- 0;
	cov_s[1:N, 1:N] <- cov_GPL2(dmat[1:N, 1:N], rhosq, etasq, 0.0000001)
	s[1:N] ~ dmnorm(mu_s[1:N], cov = cov_s[1:N, 1:N])
})

## Define Model Parmeters ----
etasq <-0.0001
rhosq <- 0.0001
r <- 0.004
curve(etasq * exp(-rhosq * x^2), from = 0, to = 100)

## Simulate Dates ----
set.seed(123)
simModel <- nimbleModel(code = spatial_demo,constants = list(start=3000,end=1800,N = N,dmat = distMat,etasq=etasq,rhosq=rhosq,r=r))
simModel$simulate('s')
simModel$simulate('local_r')
simModel$simulate('theta')
# hist(simModel$theta,xlim=c(3000,1800))
# hist(simModel$local_r)
# range(simModel$s)

## Combine Results ----
coords <- coordinates(locations) |> as.data.frame()
coords$theta  <- simModel$theta
coords$r  <- simModel$local_r
coords$s  <- simModel$s
# ggplot(coords,aes(x=x,y=y,col=s))+geom_point()+ scale_color_gradient2() + annotate('rect',xmin=0,xmax = 500000,ymin=3750000,ymax=4000000,alpha=0.1,col='black') + annotate('rect',ymin=4340000,ymax=4600000,xmin=300000,xmax=590000,alpha=0.1,col='black')
# region1 <- subset(coords,x>0&x<500000&y>3750000&y<4000000)
# region2 <- subset(coords,x>300000&x<590000&y>4340000&y<4600000)
# par(mfrow=c(2,1))
# hist(region1$theta,xlim=c(3000,1800))
# hist(region2$theta,xlim=c(3000,1800))
coords$cra_error  <-  20
coords$cra  <- round(uncalibrate(round(coords$theta))$ccCRA)
coords$medDates <- medCal(calibrate(coords$cra,coords$cra_error))

# Analyse Simulated Data ----
## Prepare Data, Inits, and Constants ----
data(intcal20)
inits=list(s=rep(N,0),etasq=0.0001,rhosq=0.001,r=0.0004,theta=coords$medDates)
inits$cov_s <- Ccov_GPL2(distMat, inits$rhosq, inits$etasq, 0.0000001)
inits$s <-  t(chol(inits$cov_s)) %*% rnorm(N)
inits$s <- inits$s[ , 1]  # so can give nimble a vector rather than one-column matrix
inits$local_r <- inits$r + inits$s
constants <- list(N=N,dists=distMat,calBP=intcal20$CalBP,C14BP=intcal20$C14Age,C14err=intcal20$C14Age.sigma,a=3000,b=1800)
data <- list(cra=coords$cra,cra_error=coords$cra_error)
## Core Model Analysis ----
runFun <- function(seed, data, constants, inits, niter, nburnin, thin) {
  library(nimbleCarbon)
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

  assign('cov_GPL2',cov_GPL2,envir=.GlobalEnv)

  model <- nimbleCode({
  for (i in 1:N){
    # Model
    local_r[i] <- r + s[i]
    theta[i] ~ dExponentialGrowth(a=a,b=b,r=local_r[i])
    
    c14age[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
    sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
    sigmaDate[i] <- (cra_error[i]^2+sigmaCurve[i]^2)^(1/2);
    cra[i] ~ dnorm(mean=c14age[i],sd=sigmaDate[i]);
  }
  #priors
  r ~ dnorm(0,sd=0.001)
  etasq ~ dexp(1);
  rhosq ~ dexp(1000);
  mu_s[1:N] <- 0;
  cov_s[1:N, 1:N] <- cov_GPL2(dists[1:N, 1:N], rhosq, etasq, 0.001)
  s[1:N] ~ dmnorm(mu_s[1:N], cov = cov_s[1:N, 1:N])
  }) 

  model <- nimbleModel(model,constants = constants,data=data,inits=inits)
  cModel <- compileNimble(model)
  conf <- configureMCMC(model)
  conf$addMonitors('s')
  conf$removeSamplers('s[1:1000]')
#  conf$addSampler(c('r','s[1:1000]'), type='AF_slice') 
  conf$addSampler(c('r','s[1:1000]'), type='RW_block') 

  MCMC <- buildMCMC(conf)
  cMCMC <- compileNimble(MCMC)
  results <- runMCMC(cMCMC, nchain=1,niter = niter, thin=thin,nburnin = nburnin,samplesAsCodaMCMC = T,setSeed = seed, progressBar = TRUE) 
  return(results)
}
# Run MCMC in Parallel----
ncores <- 3
cl <- makeCluster(ncores)
# Run the model in parallel:
seeds = c(123,456,789)
chain_output = parLapply(cl = cl, X = seeds, fun = runFun, data = data, 
constants = constants, inits = inits, niter = 100000, nburnin = 50000,thin = 5)
stopCluster(cl)
# Convert into a mcmc.list object for diagnostic (see below)
chain_output=coda::mcmc.list(chain_output)

# MCMC Diagnostics and Posterior Plots----
coda::gelman.diag(chain_output) #highest rhat is 1.01 (sigma), all random effets are 1.00




