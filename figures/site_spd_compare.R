# Load Libraries ----
library(here)
library(nimbleCarbon)
library(rcarbon)
source(here('src','unif2.R'))
load(here('data','sitedata.RData'))
load(here('data','c14data.RData'))


ricearrival = data.frame(regions=1:8,
			 m = c(-1039,-570,-910,-824,-648,-271,-152,-428),
			 hi = c(-1251,-735,-1061,-946,-754,-471,-434,-709),
			 lo = c(-872,-430,-779,-703,-560,-124,42,-203))


sitedb$RegionID = match(sitedb$Area,win.riceregion$riceregion)
sort(table(sitedb$RegionID))
nsim  <- 1000

c14db  <- subset(c14db,C14Age<= 3500 & C14Age >= 1500 & !Region%in%c('Hokkaido','Okinawa') & Material=='Terrestrial')
caldates  <- calibrate(c14db$C14Age,c14db$C14Error)
caldates.n  <- calibrate(c14db$C14Age,c14db$C14Error,normalised=FALSE)
ii  <- which.CalDates(caldates,BP<3100&BP>1650,p=0.5)
calibrated.dates  <- caldates[ii]
calibrated.dates.n <- caldates.n[ii]
c14db  <- c14db[ii,]

bin50  <- binPrep(sites = c14db$SiteID,ages=calibrated.dates,h=50)
i  <- thinDates(ages=c14db$C14Age,  errors=c14db$C14Error, bins=bin50, size=1, thresh=1,seed=123,method='splitsample')
c14db  <-  c14db[i,]
calibrated.dates  <- calibrated.dates[i]
calibrated.dates.n  <- calibrated.dates.n[i]



densMats  <- spdMats <- spdMats.n <-   vector('list',length=8)

for (k in 1:8)
{
	# Settlement Data 
	tmp  <- subset(sitedb,RegionID==k)
	a  <- 3200
	b  <- 1750
	densMat  <- matrix(NA,nrow=nsim,ncol=length(a:b))

	for (i in 1:nsim)
	{
		d  <- runif(nrow(tmp),min=tmp$end,max=tmp$start)
		dens  <- density(d,bw=50)
		densMat[i,]  <- approx(x = dens$x,y=dens$y,xout=a:b)$y
	}

	densMats[[k]]  <- data.frame(BP=a:b,m=apply(densMat,2,mean),lo=apply(densMat,2,quantile,0.025,na.rm=T),hi=apply(densMat,2,quantile,0.975,na.rm=T))
	# Radiocarbon Data
	j  <- which(c14db$RiceRegion==as.character(as.roman(k)))
	spdMats[[k]]  <- spd(calibrated.dates[j],timeRange=c(3200,1750),spdnormalised=T)
	spdMats.n[[k]]  <- spd(calibrated.dates.n[j],timeRange=c(3200,1750),spdnormalised=T)
}


par(mfrow=c(2,4))

for (i in 1:8)
{
	spd1  <- spdMats[[i]][[2]]
	spd2  <- spdMats.n[[i]][[2]]
	kdesites  <- densMats[[i]]
	maxSPD  <- max(c(spd1[,2],spd2[,2]),kdesites$hi)

	plot(NULL,xlim=c(3200,1750),ylim=c(0,maxSPD),axes=F,xlab='BC/AD',ylab='',yaxs='i')
	rect(xleft=BCADtoBP(ricearrival$hi[i]),xright=BCADtoBP(ricearrival$lo[i]),ybottom=-100,ytop=100,border=NA,col=adjustcolor('steelblue',alpha.f=0.3))
	polygon(c(3200:1750,1750:3200),c(spd2[,2],rep(0,nrow(spd2))),border=NA,col='lightgrey')
	lines(3200:1750,spd1[,2],lty=1)
	polygon(x=c(kdesites$BP,rev(kdesites$BP)),y=c(kdesites$lo,rev(kdesites$hi)),border=NA,col=adjustcolor('darkorange',alpha.f = 0.4))
	abline(v=BCADtoBP(ricearrival$m[i]),lty=4)
	prettyBCAD  <- pretty(BPtoBCAD(3200:1750))
	prettyBCAD[which(prettyBCAD==0)] <- 1
	axis(1,at=BCADtoBP(prettyBCAD),labels=abs(prettyBCAD))
	axis(1,at=BCADtoBP(c(seq(-2000,-100,100),1,seq(100,1500,100))),labels=NA,tck=-0.01)
	axis(2)
	axis(3,at=c(0,-500,-750,-1000)+BCADtoBP(ricearrival$m[i]),labels=c(0,500,750,1000),padj=+.5)
	mtext('Years from arrival of rice farming',side=3,line=2)
	box()
	legend('topleft',legend=paste0(as.roman(i)),cex=2,bty='n')
}



