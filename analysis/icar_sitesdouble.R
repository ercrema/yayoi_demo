library(ggplot2)
library(dplyr)
library(gridExtra)
library(nimbleCarbon)
library(rcarbon)
library(rnaturalearth)
library(sf)
library(maptools)
library(rgeos)
library(spdep)
library(here)
library(parallel)
# Source Uniform model ----
source(here('src','unif2.R'))

# Load and prepare 14C and Spatial Data ----
load(here('data','sitedata.RData'))
load(here('data','c14data.RData'))
timespan  <- 500
sitedb$a = sitedb$ricearrival + timespan
sitedb$b = sitedb$ricearrival - timespan

# Subset based on probability mass ----
foo  <- function(x1,y1,x2,y2)
{
	d  <- x1 - y1
	p  <- 1/d
	if (x1<=x2 & y1>=y2) {return(1)}
	if (x1>x2 & y1>x2) {return(0)}
	if (x1<y2 & y1<y2) {return(0)}
	if (x1>=x2 & y1>y2) {return(abs(x2-y1)*p)}
	if (x1<x2 & y1<y2) {return(abs(x1-y2)*p)}
	if (x1>=x2 & y1<=y2) {return(abs(x2-y2)*p)}
}

# Compute probabilities within window
sitedb$prob.within  <- unlist(apply(cbind(sitedb$start,sitedb$end,sitedb$a,sitedb$b),1,function(x){foo(x[1],x[2],x[3],x[4])}))
sitedb <- subset(sitedb,prob.within>0.5)

# Define total number of sites
N <- nrow(sitedb)

# Define Spatial Structure
# Create RiceRegion ID and add to c14db
sitedb$RegionID = match(sitedb$Area,win.riceregion$riceregion)
sort(table(sitedb$RegionID))

# Define Data
d <- list(X=round(sitedb$m),s=round(sitedb$s))

# Define Constants
constants  <-  list(Nregions=Nregions,N=N,adj=nbInfo.rice$adj,weights=nbInfo.rice$weights,num=nbInfo.rice$num,L=length(nbInfo.rice$adj))
constants$a <- sitedb$a
constants$b <- sitedb$b
constants$chp <- sitedb$ricearrival
constants$id.region <-sitedb$RegionID

# Define theta init
theta.init <- d$X


# Core Analysis ----

runFun <- function(seed, d, constants, theta.init, nburnin, niter, thin)
{
	library(nimbleCarbon)
	## ICAR  Model
	icarmodel <- nimbleCode({
		for (i in 1:N)
		{
			theta[i] ~ dDoubleExponentialGrowth(a=a[i],b=b[i],r1=s1[id.region[i]],r2=s2[id.region[i]],mu=chp[i])
			X[i] ~ dunif2(m=theta[i],s=s[i])
		}

		# ICAR Model Prior
		s1[1:Nregions] ~ dcar_normal(adj[1:L], weights[1:L], num[1:Nregions], tau1, zero_mean =0)
		s2[1:Nregions] ~ dcar_normal(adj[1:L], weights[1:L], num[1:Nregions], tau2, zero_mean =0)
		tau1 <- 1/sigma1^2
		tau2 <- 1/sigma2^2
		sigma1 ~ dunif(0,100)
		sigma2 ~ dunif(0,100)
	})

	## Setup Init
	set.seed(seed)
	inits <- list(sigma1 = runif(1,0,100),sigma2=runif(1,0,100), s1 = rnorm(constants$Nregions,sd=0.001),s2 = rnorm(constants$Nregions,sd=0.001),theta=theta.init)

	#MCMC
	model <- nimbleModel(icarmodel, constants = constants, data = d, inits = inits)
	cModel <- compileNimble(model)
	conf <- configureMCMC(model, monitors = c('s1','s2','sigma1','sigma2'))
	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC, project = cModel)
	samples <- runMCMC(cMCMC, niter = niter, thin=thin,nburnin = nburnin,samplesAsCodaMCMC = T, setSeed = seed)
	return(samples)
}

# Run MCMC ----
ncores  <-  4
cl <- makeCluster(ncores)
seeds  <-  c(12,34,56,78)
niter  <- 250000
nburnin  <- 125000
thin  <- 5

chain_output  <- parLapply(cl = cl, X = seeds, fun = runFun, d = d,constants = constants, theta = theta.init, niter = niter, nburnin = nburnin,thin = thin)
stopCluster(cl)
icar.samples <- coda::mcmc.list(chain_output)
rhats.sitesdouble  <- coda::gelman.diag(icar.samples)
icar.sitesdouble  <- do.call(rbind.data.frame,icar.samples)
save(rhats.before500,icar.before500,file=here('results','icar_sitesdoubleRes.RData'))
