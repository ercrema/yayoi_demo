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

# Load and prepare 14C and Spatial Data ----
load(here('data','c14data.RData'))
timespan  <- 500
c14db$fromFarming  <- c14db$C14Age - c14db$ricearrival
c14db$a = c14db$ricearrival + timespan
c14db$b = c14db$ricearrival - timespan
c14db = subset(c14db, C14Age < a+500 & C14Age > b-500 & Material == 'Terrestrial' & !PrefectureNameEn %in% c('Hokkaido','Okinawa')) |> select(C14Age,C14Error,a,b,ricearrival,SiteID,Prefecture=PrefectureNameEn,Longitude,Latitude,RiceRegion)

regions  <- c('I','II','III','IV','V','VI','VII','VIII')
regionList  <- vector('list',length=length(regions))
for (i in 1:length(regions))
{
	tmp.subset  <- subset(c14db,RiceRegion==regions[i])
	tmp.cal  <- calibrate(tmp.subset$C14Age,tmp.subset$C14Error)
	tmp.subset$medCal  <- medCal(tmp.cal)
	tmp.subset$include  <- FALSE
	a  <-  tmp.subset$a[1]
	b  <- tmp.subset$b[1]
	tmp.subset$include[which.CalDates(tmp.cal,BP<a&BP>b,p=0.5)]  <- TRUE
	regionList[[i]]  <- tmp.subset
}

c14db  <- do.call(rbind.data.frame,regionList) |> subset(x=_,include==TRUE)
calibrated.dates  <- calibrate(c14db$C14Age,c14db$C14Error)
# Site Level Thinning
bin50  <- binPrep(sites = c14db$SiteID,ages=calibrated.dates,h=50)
i  <- thinDates(ages=c14db$C14Age,  errors=c14db$C14Error, bins=bin50, size=1, thresh=1,seed=123,method='splitsample')
c14db  <-  c14db[i,]

# Define initialisation values for theta
c14db$theta <- c14db$medCal
c14db$theta  <- ifelse(c14db$theta>=c14db$a,c14db$a-1,c14db$theta)
c14db$theta  <- ifelse(c14db$theta<=c14db$b,c14db$b+1,c14db$theta)
# any(c14db$theta>c14db$a)
# any(c14db$theta<c14db$b)
theta.init <- c14db$theta

# Define total number of dates
N <- nrow(c14db)

# Define Spatial Structure
# Create RiceRegion ID and add to c14db
c14db$RegionID = match(c14db$RiceRegion,win.riceregion$riceregion)
sort(table(c14db$RegionID))

# Define Data
d <- list(cra=c14db$C14Age,cra_error=c14db$C14Error)

# Define Constants
data(intcal20)
constants  <-  list(Nregions=Nregions,N=N,adj=nbInfo.rice$adj,weights=nbInfo.rice$weights,num=nbInfo.rice$num,L=length(nbInfo.rice$adj),calBP=intcal20$CalBP,C14BP=intcal20$C14Age,C14err=intcal20$C14Age.sigma)
constants$a <- c14db$a
constants$b <- c14db$b
constants$chp <- c14db$ricearrival
constants$id.region <-c14db$RegionID



# Core Analysis ----

runFun <- function(seed, d, constants, theta.init, nburnin, niter, thin)
{
	library(nimbleCarbon)
	## ICAR  Model
	icarmodel <- nimbleCode({ 
		for (i in 1:N)
		{
			theta[i] ~ dDoubleExponentialGrowth(a=a[i],b=b[i],r1=s1[id.region[i]],r2=s2[id.region[i]],mu=chp[i])
			c14age[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
			sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
			sigmaDate[i] <- (cra_error[i]^2+sigmaCurve[i]^2)^(1/2);
			cra[i] ~ dnorm(mean=c14age[i],sd=sigmaDate[i]);
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
rhats.c14double  <- coda::gelman.diag(icar.samples)
icar.c14double <- do.call(rbind.data.frame,icar.samples)
save(rhats.before500,icar.before500,file=here('results','icar_c14doubleRes.RData'))
