library(here)
library(nimbleCarbon)
library(dplyr)
library(rcarbon)
library(parallel)

# Load Data ----
load(here('data','c14data.RData'))
c14db = subset(c14db,C14Age<4000 & C14Age> 300 &  Material == 'Terrestrial' & anthropic == TRUE) |> select(C14Age,C14Error,SiteID,Region=Region3,RegionID=RegionID3) |> arrange(RegionID)


# Site Level Thinning ----
calibrated.dates <- calibrate(c14db$C14Age,c14db$C14Error)
bin  <- binPrep(sites = c14db$SiteID,ages=calibrated.dates,h=50)
i  <- thinDates(ages=c14db$C14Age,  errors=c14db$C14Error, bins=bin, size=1, thresh=1,seed=123,method='splitsample')
c14db  <-  c14db[i,]
calibrated.dates <- calibrated.dates[i]

# Setup Regions ----
regionStartDates = unique(select(c14db,Region,RegionID)) |> arrange(RegionID)
regionStartDates$a1 = BCADtoBP(c(-1283,-1117,-981,-833,-685,-515,-362)) # 1283, 1117, 981, 833, 685, 515, 362
regionStartDates$a2 = BCADtoBP(c(-900))
regionStartDates$b1 = regionStartDates$a1 - 1000
regionStartDates$b2 = BCADtoBP(250)

# Subset Data to Relevant Temporal Window of Analyses ----
c14db$med = medCal(calibrated.dates)
c14db$r1 = NA
c14db$r2 = NA

for (i in 1:nrow(regionStartDates))
{
	a1 = regionStartDates$a1[i]
	a2 = regionStartDates$a2[i]
	b1 = regionStartDates$b1[i]
	b2 = regionStartDates$b2[i]
	ii = which(c14db$RegionID == i)
	c14db$r1[ii] = FALSE
	c14db$r2[ii] = FALSE
	regionalDates = calibrated.dates[ii]
	iii = ii[which.CalDates(regionalDates,BP < a1 & BP > b1,p=0.5)]
	c14db$r1[iii] = TRUE
	iii = ii[which.CalDates(regionalDates,BP < a2 & BP > b2,p=0.5)]
	c14db$r2[iii] = TRUE
}

# Define Data & Constants ----
data(intcal20)
c14db1 = subset(c14db,r1==TRUE)
c14db2 = subset(c14db,r2==TRUE)

#Constants
const1 <- const2 <- list(calBP=intcal20$CalBP,C14BP=intcal20$C14Age,C14err=intcal20$C14Age.sigma)
const1$Nregions <- const2$Nregions <- max(c14db$RegionID)
const1$Ndates <- nrow(c14db1)
const2$Ndates <- nrow(c14db2)
const1$region.id <- c14db1$RegionID
const2$region.id <- c14db2$RegionID
const1$starts <- regionStartDates$a1
const2$starts <- regionStartDates$a2
const1$ends <- regionStartDates$b1
const2$ends <- regionStartDates$b2

# Data
d1 <- list(cra=c14db1$C14Age,cra_error=c14db1$C14Error)
d2 <- list(cra=c14db2$C14Age,cra_error=c14db2$C14Error)

# Define Sensible initialisation values for theta ----
c14db1 <- left_join(c14db1,regionStartDates)
c14db2 <- left_join(c14db2,regionStartDates)
c14db1$theta.init <- c14db1$med
c14db2$theta.init <- c14db2$med
c14db1$theta.init[which(c14db1$theta.init>=c14db1$a1)]  <- c14db1$a1[which(c14db1$theta.init>=c14db1$a1)] - 1
c14db1$theta.init[which(c14db1$theta.init<=c14db1$b1)] <- c14db1$b1[which(c14db1$theta.init<=c14db1$b1)] + 1
c14db2$theta.init[which(c14db2$theta.init>=c14db2$a)]  <- c14db2$a2[which(c14db2$theta.init>=c14db2$a2)] - 1
c14db2$theta.init[which(c14db2$theta.init<=c14db2$b2)] <- c14db2$b2[which(c14db2$theta.init<=c14db2$b2)] + 1
theta.init.1 <- c14db1$theta.init
theta.init.2  <- c14db2$theta.init

# Sanity Check Plot
 spd1 <- spd2 <- vector('list',length=nrow(regionStartDates))
 for (i in 1:nrow(regionStartDates))
 {
 	ii <- which(c14db1$RegionID==i)
 	cal.1 <- calibrate(c14db1$C14Age[ii],c14db1$C14Error[ii],normalised=T)
 	spd1[[i]] <- spd(cal.1,timeRange=c(regionStartDates$a1[i],regionStartDates$b1[i]))
 
 	ii <- which(c14db2$RegionID==i)
 	cal.2 <- calibrate(c14db2$C14Age[ii],c14db2$C14Error[ii],normalised=T)
 	spd2[[i]] <- spd(cal.2,timeRange=c(regionStartDates$a2[i],regionStartDates$b2[i]))
 }
 
 par(mfrow=c(7,1),mar=c(4,1,1,1))
 lapply(spd2,plot)
 lapply(spd1,plot)



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

chain_output_1000yrs = parLapply(cl = cl, X = seeds, fun = runFun, d = d1,constants = const1, theta = theta.init.1, niter = niter, nburnin = nburnin,thin = thin)
chain_output_yayoi = parLapply(cl = cl, X = seeds, fun = runFun, d = d2,constants = const2, theta = theta.init.2, niter = niter, nburnin = nburnin,thin = thin)
stopCluster(cl)
post.sample.1000yrs=coda::mcmc.list(chain_output_1000yrs)
post.sample.yayoi=coda::mcmc.list(chain_output_yayoi)
coda::gelman.diag(post.sample.1000yrs)
coda::gelman.diag(post.sample.yayoi)

save(post.sample.1000yrs,file=here('R_images','regdemo_1000yrs.RData'))
save(post.sample.yayoi,file=here('R_images','regdemo_yayoi.RData'))
