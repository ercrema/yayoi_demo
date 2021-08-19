library(here)
library(nimbleCarbon)
library(dplyr)
library(rcarbon)
library(parallel)

# Load Data ----
load(here('data','c14data.RData'))
c14db = subset(c14db,C14Age<4000 & C14Age> 300 &  Material == 'Terrestrial') |> select(C14Age,C14Error,SiteID,Region=Region3,RegionID=RegionID3) |> arrange(RegionID)


# Site Level Thinning ----
calibrated.dates <- calibrate(c14db$C14Age,c14db$C14Error)
bin50  <- binPrep(sites = c14db$SiteID,ages=calibrated.dates,h=50)
i  <- thinDates(ages=c14db$C14Age,  errors=c14db$C14Error, bins=bin50, size=1, thresh=1,seed=123,method='splitsample')
c14db  <-  c14db[i,]
calibrated.dates <- calibrated.dates[i]

# Setup Regions ----
regionStartDates = unique(select(c14db,Region,RegionID)) |> arrange(RegionID)
# regionStartDates$a = 3000
regionStartDates$a = BCADtoBP(c(-1283,-1117,-981,-833,-685,-515,-362)) # 1283, 1117, 981, 833, 685, 515, 362
regionStartDates$b500 = regionStartDates$a - 500
regionStartDates$b1000 = regionStartDates$a - 1000

# Subset Data to Relevant Temporal Window of Analyses ----
c14db$med = medCal(calibrated.dates)
c14db$r500 = NA
c14db$r1000 = NA

for (i in 1:nrow(regionStartDates))
{
	a = regionStartDates$a[i]
	b500 = regionStartDates$b500[i]
	b1000 = regionStartDates$b1000[i]
	ii = which(c14db$RegionID == i)
	c14db$r500[ii] = FALSE
	c14db$r1000[ii] = FALSE
	regionalDates = calibrated.dates[ii]
	iii = ii[which.CalDates(regionalDates,BP < a & BP > b500,p=0.5)]
	c14db$r500[iii] = TRUE
	iii = ii[which.CalDates(regionalDates,BP < a & BP > b1000,p=0.5)]
	c14db$r1000[iii] = TRUE
}

# Define Data & Constants ----
data(intcal20)
c14db500 = subset(c14db,r500==TRUE)
c14db1000 = subset(c14db,r1000==TRUE)

#Constants
const500 <- const1000 <- list(calBP=intcal20$CalBP,C14BP=intcal20$C14Age,C14err=intcal20$C14Age.sigma)
const500$Nregions <- const1000$Nregions <- max(c14db$RegionID)
const500$Ndates <- nrow(c14db500)
const1000$Ndates <- nrow(c14db1000)
const500$region.id <- c14db500$RegionID
const1000$region.id <- c14db1000$RegionID
const500$starts <- const1000$starts <- regionStartDates$a
const500$ends <- regionStartDates$b500
const1000$ends <- regionStartDates$b1000

# Data
d500 <- list(cra=c14db500$C14Age,cra_error=c14db500$C14Error)
d1000 <- list(cra=c14db1000$C14Age,cra_error=c14db1000$C14Error)

# Define Sensible initialisation values for theta ----
c14db500 <- left_join(c14db500,regionStartDates)
c14db1000 <- left_join(c14db1000,regionStartDates)
c14db500$theta.init <- c14db500$med
c14db1000$theta.init <- c14db1000$med
c14db500$theta.init[which(c14db500$theta.init>=c14db500$a)]  <- c14db500$a[which(c14db500$theta.init>=c14db500$a)] - 1
c14db500$theta.init[which(c14db500$theta.init<=c14db500$b500)] <- c14db500$b[which(c14db500$theta.init<=c14db500$b500)] + 1
c14db1000$theta.init[which(c14db1000$theta.init>=c14db1000$a)]  <- c14db1000$a[which(c14db1000$theta.init>=c14db1000$a)] - 1
c14db1000$theta.init[which(c14db1000$theta.init<=c14db1000$b1000)] <- c14db1000$b[which(c14db1000$theta.init<=c14db1000$b1000)] + 1
theta.init.500 <- c14db500$theta.init
theta.init.1000  <- c14db1000$theta.init

# Sanity Check Plot
 spd.500 <- spd.1000 <- vector('list',length=nrow(regionStartDates))
 for (i in 1:nrow(regionStartDates))
 {
 	ii <- which(c14db500$RegionID==i)
 	cal.500 <- calibrate(c14db500$C14Age[ii],c14db500$C14Error[ii],normalised=F)
 	spd.500[[i]] <- spd(cal.500,timeRange=c(regionStartDates$a[i],regionStartDates$b500[i]))
 
 	ii <- which(c14db1000$RegionID==i)
 	cal.1000 <- calibrate(c14db1000$C14Age[ii],c14db1000$C14Error[ii],normalised=F)
 	spd.1000[[i]] <- spd(cal.1000,timeRange=c(regionStartDates$a[i],regionStartDates$b1000[i]))
 }
 
 par(mfrow=c(3,3))
 lapply(spd.1000,plot)



# Core Bayesian Analysis ----

runFun <- function(seed, d, constants, theta.init, nburnin, niter, thin)
{
	# Load Package
	library(nimbleCarbon)
	# Define core model
	growthmodel <- nimbleCode({
		for (i in 1:Ndates)
		{
			theta[i] ~ dExponentialGrowth(a=starts[region.id[i]],b=ends[region.id[i]],r=r[region.id[i]]);
			mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
			sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
			error[i] <- (cra_error[i]^2+sigmaCurve[i]^2)^(1/2);
			cra[i] ~ dnorm(mean=mu[i],sd=error[i]);
		}

		# prior
		for (j in 1:Nregions)
		{
			r[j] ~ dnorm(rbar,rsd)
		}

		# hyperpriors
		rbar ~ dnorm(0.0004,0.05)
		rsd ~ dexp(1/0.001)
	})

	# Set Inits
	set.seed(seed)
	inits <- list(theta=theta.init)
	inits$rbar  <- rnorm(1,0.0004,0.05)
	inits$rsd  <- rexp(1,1/0.001)
	inits$r  <- rnorm(constants$Nregions,inits$rbar,inits$rsd)

	# Run MCMC
	model <- nimbleModel(growthmodel,constants = constants,data=d,inits=inits)
	cModel <- compileNimble(model)
	conf <- configureMCMC(model)
	conf$addMonitors('r')
	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC)
	results <- runMCMC(cMCMC, niter = niter, thin=thin,nburnin = nburnin,samplesAsCodaMCMC = T, setSeed = seed)
	return(results)
}


# Run MCMC ----
ncores  <-  3
cl <- makeCluster(ncores)
seeds  <-  c(123,456,789)
niter  <- 20000
nburnin  <- 10000
thin  <- 1

chain_output_500 = parLapply(cl = cl, X = seeds, fun = runFun, d = d500,constants = const500, theta = theta.init.500, niter = niter, nburnin = nburnin,thin = thin)
chain_output_1000 = parLapply(cl = cl, X = seeds, fun = runFun, d = d1000,constants = const1000, theta = theta.init.1000, niter = niter, nburnin = nburnin,thin = thin)
stopCluster(cl)
post.sample500=coda::mcmc.list(chain_output_500)
post.sample1000=coda::mcmc.list(chain_output_1000)
coda::gelman.diag(post.sample500)
coda::gelman.diag(post.sample1000)

save(post.sample500,file=here('R_images','regdemo500.RData'))
save(post.sample1000,file=here('R_images','regdemo1000.RData'))
