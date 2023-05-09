# Load libraries and R images ----
library(here)
library(nimbleCarbon)
library(rcarbon)

load(here('data','c14data.RData'))
load(here('results','icar_c14doubleRes500.RData'))

# Thinning
calibrated.dates  <- calibrate(c14db$C14Age,c14db$C14Error)
bin50  <- binPrep(sites = c14db$SiteID,ages=calibrated.dates,h=50)
i  <- thinDates(ages=c14db$C14Age,  errors=c14db$C14Error, bins=bin50, size=1, thresh=1,seed=123,method='splitsample')
c14db  <-  c14db[i,]

# Sim settings
nsim  <- 500
ncores  <- 5
ppcheck500  <- vector("list",length=8)

# PPcheck ----
for (i in 1:8)
{
	print(paste('Running region',i))
	tmp.c14data  <- subset(c14db,RiceRegion==as.character(as.roman(i)))
	params  <- list()
	ra  <- unique(tmp.c14data$ricearrival)
	a  <-  ra + 500
	b  <- ra - 500
	params$r1  <- icar.c14double500[,i] 
	params$r2  <- icar.c14double500[,i+8]
        params$mu  <- rep(ra,nrow(icar.c14double500))
	tmp.res  <- postPredSPD(tmp.c14data$C14Age,tmp.c14data$C14Error,calCurve='intcal20',model=dDoubleExponentialGrowth,a=a,b=b,params=params,nsim=nsim,ncores=ncores,verbose=FALSE,method='uncalsample')
	ppcheck500[[i]] <- tmp.res
	assign(paste0('ppcheck.500_',i),tmp.res)
}

save(ppcheck500,file=here('results','ppcheck500.RData'))

save(ppcheck.500_1,file=here('results','ppcheck500_1.RData'))
save(ppcheck.500_2,file=here('results','ppcheck500_2.RData'))
save(ppcheck.500_3,file=here('results','ppcheck500_3.RData'))
save(ppcheck.500_4,file=here('results','ppcheck500_4.RData'))
save(ppcheck.500_5,file=here('results','ppcheck500_5.RData'))
save(ppcheck.500_6,file=here('results','ppcheck500_6.RData'))
save(ppcheck.500_7,file=here('results','ppcheck500_7.RData'))
save(ppcheck.500_8,file=here('results','ppcheck500_8.RData'))
