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
	assign(paste0('ppcheck.750_',i),tmp.res)
}

save(ppcheck.750_1,file=here('results','ppcheck750_1.RData'))
save(ppcheck.750_2,file=here('results','ppcheck750_2.RData'))
save(ppcheck.750_3,file=here('results','ppcheck750_3.RData'))
save(ppcheck.750_4,file=here('results','ppcheck750_4.RData'))
save(ppcheck.750_5,file=here('results','ppcheck750_5.RData'))
save(ppcheck.750_6,file=here('results','ppcheck750_6.RData'))
save(ppcheck.750_7,file=here('results','ppcheck750_7.RData'))
save(ppcheck.750_8,file=here('results','ppcheck750_8.RData'))
