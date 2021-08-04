library(ggplot2)
library(nimbleCarbon)
library(rcarbon)
library(rnaturalearth)
library(maptools)
library(rgeos)
library(spdep)
library(here)

# Setup Window of Analysis and Site Locations ----

## Generate Window
win <- ne_states(country = "japan") |> subset(!name_vi %in% c("Okinawa","Hokkaido")) |> spTransform(CRSobj = CRS("+proj=utm +zone=54 ellps=WGS84"))
## Create Prefecture ID
Npref <- length(win)
win@data$ID = 1:Npref

## Compute Distance Matrices
pref.centers  <-  gCentroid(win,byid = TRUE)
distMatrix  <- coordinates(pref.centers) |> dist() |> as.matrix()
distMatrix  <- distMatrix/1000



## Define Neighbourhood Structure
W_nb <- poly2nb(win, row.names =  win@data$ID)
## Determine neighborhood/adjacency information needed for neighborhood-based CAR model
nbInfo <- nb2WB(W_nb)
# A vector of indices indicating which regions are neighbors of which.
# head(nbInfo$adj, n = 30)
# A vector of weights. In this case, all weights are 1.
# head(nbInfo$weights)
# A vector of length N indicating how many neighbors each region has.
# This helps map the adj vector to each region.
# nbInfo$num


# Simulate Data (using GP) ----
## Generate Samples Sizes for Each Region
avgSample = 1000
nSamples <- rpois(Npref,lambda=avgSample/Npref)
N = sum(nSamples)
id.pref  <- rep(1:Npref,nSamples)

## USE GP to generate observed data
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
	for (i in 1:Npref){
		local_r[i]  <- r + s[i]
	}
	for (i in 1:N){
		theta[i] ~ dExponentialGrowth(a=start,b=end,r=local_r[id.pref[i]])
	}
	mu_s[1:Npref] <- 0;
	cov_s[1:Npref, 1:Npref] <- cov_GPL2(dmat[1:Npref, 1:Npref], rhosq, etasq, 0.0000001)
	s[1:Npref] ~ dmnorm(mu_s[1:Npref], cov = cov_s[1:Npref, 1:Npref])
})

### Parameters
etasq <-0.00001
rhosq <- 0.00006
r <- 0.004
curve(etasq * exp(-rhosq * x^2), from = 0, to = 200)


### Actual Simulation
set.seed(1234)
simModel <- nimbleModel(code = spatial_demo,constants = list(start=3000,end=1800,N = N,dmat = distMatrix,etasq=etasq,rhosq=rhosq,r=r,Npref=Npref,id.pref=id.pref))
simModel$simulate('s')
simModel$simulate('local_r')
simModel$simulate('theta')

win@data$r = simModel$local_r
win.sf  <-  as(win,'sf')
win.sf  <- cbind(win.sf,st_coordinates(st_centroid(win.sf)))
ggplot(win.sf,aes(fill=r))  + geom_sf() + scale_fill_gradient2()


d  <- data.frame(ID=1:N,id.pref=id.pref,cra=round(uncalibrate(round(simModel$theta))$ccCRA),cra_error=20)
d$medDates  <- medCal(calibrate(d$cra,d$cra_error))
d$thetainit  <- d$medDates
a <- 3000
b <- 1800
d$thetainit[which(d$thetainit>=a)]  <- a - 1
d$thetainit[which(d$thetainit<=b)]  <- b + 1



# Analyse Data ----
## Core Simulation Model
icarmodel <- nimbleCode({
	for (i in 1:N)
	{
		theta[i] ~ dExponentialGrowth(a=a,b=b,r=s[id.pref[i]])
		c14age[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
		sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
		sigmaDate[i] <- (cra_error[i]^2+sigmaCurve[i]^2)^(1/2);
		cra[i] ~ dnorm(mean=c14age[i],sd=sigmaDate[i]);
	}

	# ICAR Model Prior
	s[1:Npref] ~ dcar_normal(adj[1:L], weights[1:L], num[1:Npref], tau, zero_mean =0)
	tau <- 1/sigma^2
	sigma ~ dunif(0,100)
})

## Setup Constants
data(intcal20)
constants  <-  list(a=a,b=b,Npref=Npref,N=N,adj=nbInfo$adj,weights=nbInfo$weights,num=nbInfo$num,L=length(nbInfo$adj),id.pref=id.pref,calBP=intcal20$CalBP,C14BP=intcal20$C14Age,C14err=intcal20$C14Age.sigma)

## Setup Data
data <- list(cra=d$cra,cra_error=d$cra_error) 

## Setup Init
inits <- list(sigma = 0.001, s = rnorm(Npref,sd=0.001),theta=d$thetainit)


## Run
model <- nimbleModel(icarmodel, constants = constants, data = data, inits = inits)
cModel <- compileNimble(model)
conf <- configureMCMC(model, monitors = c('s','sigma'))
conf$printSamplers()
MCMC <- buildMCMC(conf)
cMCMC <- compileNimble(MCMC, project = cModel)
samples <- runMCMC(cMCMC, niter = 10000, nburnin = 5000, thin = 1)
samples.mu <- apply(samples,2,median)
plot(simModel$local_r,samples.mu[1:45])
abline(a=0,b=1)
