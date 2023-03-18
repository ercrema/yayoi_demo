library(here)  
library(nimbleCarbon)
library(rcarbon)
library(dplyr)
library(parallel)

# Load Results and Raw Data ----
# Load Results MCMC
# load(here('R_images','regdemo_1000yrs.RData'))
load(here('R_images','regdemo_800yrs.RData'))
load(here('R_images','regdemo_yayoi.RData'))


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
gap  <- 800
regionStartDates = unique(select(c14db,Region,RegionID)) |> arrange(RegionID)
regionStartDates$a1 = BCADtoBP(c(-1283,-1117,-981,-833,-685,-515,-362)) # 1283, 1117, 981, 833, 685, 515, 362
regionStartDates$a2 = BCADtoBP(c(-900))
regionStartDates$b1 = regionStartDates$a1 - gap
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


# Posterior Predictive Check ----
ppcheck_yayoi <- ppcheck_800yrs <- vector('list',length=nrow(regionStartDates))

for (i in 1:nrow(regionStartDates)) 
{
	print(paste0('Running posterior check; region ',i,' of ',nrow(regionStartDates)))
	obs.cra <- d1$cra[c14db1$RegionID==i]
	obs.error <- d1$cra_error[c14db1$RegionID==i]
	a <- regionStartDates$a1[i]
	b <- regionStartDates$b1[i]
	params <- list(r=post.sample.800yrs[[1]][,paste0('r[',i,']')])
	ppcheck_800yrs[[i]] <- postPredSPD(obs.cra,errors=obs.error,calCurve = 'intcal20',model = dExponentialGrowth,a = a,b=b,params=params,method='uncalsample',nsim = 500,ncores = 2,verbose=TRUE)

	# obs.cra <- d2$cra[c14db2$RegionID==i]
	# obs.error <- d2$cra_error[c14db2$RegionID==i]
	# a <- regionStartDates$a2[i]
	# b <- regionStartDates$b2[i]
	# params <- list(r=post.sample.yayoi[[1]][,paste0('r[',i,']')])
	# ppcheck_yayoi[[i]] <- postPredSPD(obs.cra,errors=obs.error,calCurve = 'intcal20',model = dExponentialGrowth,a = a,b=b,params=params,method='uncalsample',nsim = 500,ncores = 3,verbose=TRUE)
}


save(ppcheck_yayoi,file=here('R_images','ppcheck_yayoi.RData'))
save(ppcheck_800yrs,file=here('R_images','ppcheck_800yrs.RData'))


