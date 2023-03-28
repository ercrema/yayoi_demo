library(rcarbon)
library(here)
library(dplyr)
load(here('data','c14data.RData'))

a <- 7000 
b <- 1700 


# Load and prepare 14C and Spatial Data ----
# Load Data
load(here('data','c14data.RData'))
c14db$MacroRegion  <- 'SW Japan'
c14db$MacroRegion[which(c14db$Region%in%c('Kanto','Chubu'))]  <- 'Central Japan'
c14db$MacroRegion[which(c14db$Region=='Tohoku')]  <- 'NE Japan'
c14db = subset(c14db,C14Age< (a+500) & C14Age> (b-500) &  Material == 'Terrestrial' & !PrefectureNameEn %in% c('Hokkaido','Okinawa')) |> select(C14Age,C14Error,SiteID,MacroRegion,RiceRegion)

calibrated.dates  <- calibrate(c14db$C14Age,c14db$C14Error)

bin50  <- binPrep(sites = c14db$SiteID,ages=calibrated.dates,h=50)
i  <- thinDates(ages=c14db$C14Age,  errors=c14db$C14Error, bins=bin50, size=1, thresh=1,seed=123,method='splitsample')
c14db  <-  c14db[i,]
calibrated.dates <- calibrated.dates[i]

sw.japan.dates  <- calibrated.dates[which(c14db$MacroRegion=='SW Japan')]
c.japan.dates  <- calibrated.dates[which(c14db$MacroRegion=='Central Japan')]
ne.japan.dates  <- calibrated.dates[which(c14db$MacroRegion=='NE Japan')]

sw.japan.spd  <- spd(sw.japan.dates,timeRange=c(7000,1700),runm=50)
c.japan.spd  <- spd(c.japan.dates,timeRange=c(7000,1700),runm=50)
ne.japan.spd  <- spd(ne.japan.dates,timeRange=c(7000,1700),runm=50)

pdf(here('figures','spd_compare.pdf'),width=8,height=10)
par(mfrow=c(3,1),mar=c(5.2,4,1.2,1))
plot(ne.japan.spd,calendar='BCAD')
abline(v=-428,lty=5)
barCodes(x=BPtoBCAD(medCal(ne.japan.dates)),width = 10,yrng = c(0.95*par('usr')[4],par('usr')[4]))
legend('bottomright',legend=c(paste0('n=',length(ne.japan.dates))),bty='n',cex=1.5)
legend('left',legend=c('NE Japan'),cex=2,bty='n')

plot(c.japan.spd,calendar='BCAD')
abline(v=-648,lty=5)
barCodes(x=BPtoBCAD(medCal(c.japan.dates)),width = 10,yrng = c(0.95*par('usr')[4],par('usr')[4]))
legend('bottomright',legend=c(paste0('n=',length(c.japan.dates))),bty='n',cex=1.5)
legend('left',legend=c('Central Japan'),cex=2,bty='n')

plot(sw.japan.spd,calendar='BCAD')
abline(v=-1039,lty=5)
barCodes(x=BPtoBCAD(medCal(sw.japan.dates)),width = 10,yrng = c(0.95*par('usr')[4],par('usr')[4]))
legend('bottomright',legend=c(paste0('n=',length(sw.japan.dates))),bty='n',cex=1.5)
legend('left',legend=c('SW Japan'),cex=2,bty='n')
dev.off()


# Arrival Region plots ----
arrivalDF = data.frame(lo=c(-1251,-735,-1061,-946,-754,-471,-434,-709),
		       mi=c(-1039,-570,-910,-824,-648,-271,-152,-428),
		       hi=c(-872,-430,-779,-703,-560,-124,42,-203))



regionNames  <- c('I','II','III','IV','V','VI','VII','VIII')
regionCKDE  <- vector('list',length=length(regionNames))
for (i in 1:length(regionNames))
{
	tmp  <- calibrated.dates[which(c14db$RiceRegion==regionNames[i])]
	s.tmp  <- sampleDates(tmp,nsim=1000,boot=TRUE)
	regionCKDE[[i]]  <- ckde(s.tmp,timeRange=c(3500,1700),bw=50,normalised=TRUE)
}

pdf(here('figures','ckde_riceregion.pdf'),width=8,height=10)
par(mfrow=c(4,2),mar=c(5,4,2,1))
for (i in rev(1:length(regionNames)))
{
	print(i)
	plot(regionCKDE[[i]],type='envelope',calendar='BCAD',main=paste0('Area ',regionNames[[i]]))
# 	abline(v=arrivalDF$lo[i],lty=2)
# 	abline(v=arrivalDF$hi[i],lty=2)
# 	abline(v=arrivalDF$mi[i],lty=1)
	rect(ybottom=-100,ytop=100,xleft=arrivalDF$lo[i],xright=arrivalDF$hi[i],col=rgb(0.68,0.85,0.9,0.3),border=NA)
	abline(v=arrivalDF$mi[i],lty=1)
}
dev.off()

