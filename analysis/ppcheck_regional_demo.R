library(here)
library(nimbleCarbon)
library(rcarbon)
library(dplyr)
library(parallel)

# Load Results and Raw Data ----
# Load Results MCMC
load(here('R_images','regdemo500.RData'))
load(here('R_images','regdemo1000.RData'))

# Load Observed Regional Data (same as analysis/regional.demo.R)
load(here('data','c14data.RData'))
c14db = subset(c14db,C14Age<4000 & C14Age> 300 &  Material == 'Terrestrial') |> select(C14Age,C14Error,SiteID,Region=Region3,RegionID=RegionID3) |> arrange(RegionID)

## Site Level Thinning 
calibrated.dates <- calibrate(c14db$C14Age,c14db$C14Error)
bin100  <- binPrep(sites = c14db$SiteID,ages=calibrated.dates,h=100)
i  <- thinDates(ages=c14db$C14Age,  errors=c14db$C14Error, bins=bin100, size=1, thresh=1,seed=123,method='splitsample')
c14db  <-  c14db[i,]
calibrated.dates <- calibrated.dates[i]

# Setup Regions ----
regionStartDates = unique(select(c14db,Region,RegionID)) |> arrange(RegionID)
regionStartDates$a = BCADtoBP(c(-990,-890,-790,-650,-500,-480,-250))
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

data(intcal20)
c14db500 = subset(c14db,r500==TRUE)
c14db1000 = subset(c14db,r1000==TRUE)

# Data
d500 <- list(cra=c14db500$C14Age,cra_error=c14db500$C14Error)
d1000 <- list(cra=c14db1000$C14Age,cra_error=c14db1000$C14Error)

# Posterior Predictive Check ----
ppcheck500 <- ppcheck1000 <- vector('list',length=nrow(regionStartDates))

for (i in 1:nrow(regionStartDates))
{
	print(paste0('Running posterior check; region ',i,' of ',nrow(regionStartDates)))
	obs.cra <- d500$cra[c14db500$RegionID==i]
	obs.error <- d500$cra_error[c14db500$RegionID==i]
	a <- regionStartDates$a[i]
	b <- regionStartDates$b500[i]
	params <- list(r=post.sample500[[1]][,paste0('r[',i,']')])
	ppcheck500[[i]] <- postPredSPD(obs.cra,errors=obs.error,calCurve = 'intcal20',model = dExponentialGrowth,a = a,b=b,params=params,method='uncalsample',nsim = 100,ncores = 3,verbose=TRUE)

	obs.cra <- d1000$cra[c14db1000$RegionID==i]
	obs.error <- d1000$cra_error[c14db1000$RegionID==i]
	a <- regionStartDates$a[i]
	b <- regionStartDates$b1000[i]
	params <- list(r=post.sample1000[[1]][,paste0('r[',i,']')])
	ppcheck1000[[i]] <- postPredSPD(obs.cra,errors=obs.error,calCurve = 'intcal20',model = dExponentialGrowth,a = a,b=b,params=params,method='uncalsample',nsim = 100,ncores = 3,verbose=TRUE)
}






