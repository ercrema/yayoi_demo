# Load Libraries ----
library(here)
library(nimbleCarbon)
library(rcarbon)
load(here('data','sitedata.RData'))
load(here('data','c14data.RData'))


ricearrival = data.frame(regions=1:8,
			 m = c(-1039,-570,-910,-824,-648,-271,-152,-428),
			 hi = c(-1251,-735,-1061,-946,-754,-471,-434,-709),
			 lo = c(-872,-430,-779,-703,-560,-124,42,-203))


sitedb$RegionID = match(sitedb$Area,win.riceregion$riceregion)

# Window of display
c14db <- subset(c14db,Yayoi==TRUE)
a  <- 2950
b  <- 1650

calibrated.dates  <- calibrate(c14db$C14Age,c14db$C14Error)
bin50  <- binPrep(sites = c14db$SiteID,ages=calibrated.dates,h=50)
i  <- thinDates(ages=c14db$C14Age,  errors=c14db$C14Error, bins=bin50, size=1, thresh=1,seed=123,method='splitsample')
c14db  <-  c14db[i,]
calibrated.dates  <- calibrated.dates[i]
# calibrated.dates.n  <- calibrated.dates.n[i]


# Number of CKDE simulations
nsim  <- 1000
# Generate lists of storing results
densMats  <- spdMats <- spdMats.n <-   vector('list',length=8)

for (k in 1:8)
{
	# Settlement Data 
	tmp  <- subset(sitedb,RegionID==k)
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
	spdMats[[k]]  <- spd(calibrated.dates[j],timeRange=c(a,b),spdnormalised=T,runm=50)
# 	spdMats.n[[k]]  <- spd(calibrated.dates.n[j],timeRange=c(a,b),spdnormalised=T,runm=50)
}



pdf(file=here('figures','figure2.pdf'),width=5.3,height=2.5,pointsize=2)
cex=0.9
par(mfrow=c(2,4),mar=c(5,3.5,3,3.5))

for (i in 1:8)
{
	spd1  <- spdMats[[i]][[2]]
# 	spd2  <- spdMats.n[[i]][[2]]
	kdesites  <- densMats[[i]]
# 	maxSPD  <- max(c(spd1[,2],spd2[,2]))
	maxSPD  <- max(c(spd1[,2]))

	plot(NULL,xlim=c(2900,1700),ylim=c(0,maxSPD),axes=F,xlab='BC/AD',ylab='',yaxs='i',xaxs='i',cex.main=cex,cex.axis=cex,cex.lab=cex,cex=cex)
	rect(xleft=BCADtoBP(ricearrival$hi[i]),xright=BCADtoBP(ricearrival$lo[i]),ybottom=-100,ytop=100,border=NA,col=adjustcolor('steelblue',alpha.f=0.3),cex=cex)
	polygon(c(a:b,b:a),c(spd1[,2],rep(0,nrow(spd1))),border=1,col='lightgrey',cex=cex)
# 	lines(a:b,spd2[,2],lty=1,cex=cex)
	par(new=T)
	plot(NULL,xlim=c(2900,1700),ylim=c(0,max(kdesites[,-1])),axes=F,xlab='',ylab='',yaxs='i',xaxs='i',cex.main=cex,cex.axis=cex,cex.lab=cex,cex=cex)
	polygon(x=c(kdesites$BP,rev(kdesites$BP)),y=c(kdesites$lo,rev(kdesites$hi)),border=NA,col=adjustcolor('darkorange',alpha.f = 0.4),cex=cex)
	abline(v=BCADtoBP(ricearrival$m[i]),lty=4,cex=cex)
	prettyBCAD  <- pretty(BPtoBCAD(a:b))
	prettyBCAD[which(prettyBCAD==0)] <- 1
	axis(1,at=BCADtoBP(prettyBCAD),labels=abs(prettyBCAD),cex.axis=cex,cex.lab=cex,cex=cex)
	axis(1,at=BCADtoBP(c(seq(-2000,-100,100),1,seq(100,1500,100))),labels=NA,tck=-0.01,cex.axis=cex,cex.lab=cex,cex=cex)
	axis(2,cex.axis=cex,cex.lab=cex,cex=cex)
	axis(3,at=c(0,-250,-500,-750,-1000)+BCADtoBP(ricearrival$m[i]),labels=c(0,250,500,750,1000),padj=+.5,cex.axis=cex,cex.lab=cex)
	axis(4,cex.axis=cex,cex.lab=cex,cex=cex)
	mtext('Years from arrival of rice farming',side=3,line=2,cex.lab=cex*0.7,cex=cex*0.7)
	mtext('Summed Probability',side=2,line=2,cex.lab=cex*0.7,cex=cex*0.7)
	mtext('Density',side=4,line=2,cex.lab=cex*0.7,cex=cex*0.7)
	box()
	legend('bottomright',legend=paste0(as.roman(i)),cex=cex*2,bty='n')
}
dev.off()
