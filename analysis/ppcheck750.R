# Load libraries and R images ----
library(here)
library(nimbleCarbon)
library(rcarbon)

load(here('data','c14data.RData'))
load(here('results','icar_c14doubleRes750.RData'))

# Thinning
calibrated.dates  <- calibrate(c14db$C14Age,c14db$C14Error)
bin50  <- binPrep(sites = c14db$SiteID,ages=calibrated.dates,h=50)
i  <- thinDates(ages=c14db$C14Age,  errors=c14db$C14Error, bins=bin50, size=1, thresh=1,seed=123,method='splitsample')
c14db  <-  c14db[i,]

# Sim settings
nsim  <- 500
ncores  <- 5
ppcheck750  <- vector("list",length=8)

# PPcheck ----
for (i in 1:8)
{
	print(paste('Running region',i))
	tmp.c14data  <- subset(c14db,RiceRegion==as.character(as.roman(i)))
	params  <- list()
	ra  <- unique(tmp.c14data$ricearrival)
	a  <-  ra + 750
	b  <- ra - 750
	params$r1  <- icar.c14double750[,i] 
	params$r2  <- icar.c14double750[,i+8]
        params$mu  <- rep(ra,nrow(icar.c14double750))
	tmp.res  <- postPredSPD(tmp.c14data$C14Age,tmp.c14data$C14Error,calCurve='intcal20',model=dDoubleExponentialGrowth,a=a,b=b,params=params,nsim=nsim,ncores=ncores,verbose=FALSE,method='uncalsample')
	ppcheck750[[i]] <- tmp.res
}
save(ppcheck750,file=here('results','ppcheck750.RData'))
