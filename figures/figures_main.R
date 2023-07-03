# Load Libraries ----
library(here)
library(grDevices)
library(RColorBrewer)
library(rnaturalearth)
library(sf)
library(dplyr)
library(latex2exp)
library(nimbleCarbon)
library(rcarbon)

# Figure 1 ----
load(here('data','c14data.RData'))
win  <- ne_countries(scale=10,returnclass='sf')
japan <- ne_states(country = "japan",returnclass='sf') |> subset(!name_vi %in% c("Okinawa","Hokkaidō"))
df.pref.reg <- read.csv(here('data','prefecture_data.csv'))
japan  <- left_join(japan,df.pref.reg,by=c('name_ja'='JpNames'))
japan  <- group_by(japan,by=Area) |> summarise()
win  <- st_combine(win)

c14sites  <- subset(c14db,Yayoi==T|D500==T|D750==T) |> select(Longitude,Latitude) |> unique() |> st_as_sf(x=_,coords=c('Longitude','Latitude'),crs=4326)


cols  <- c(brewer.pal('Set2',n=8))


pdf(file=here('figures','figure1.pdf'),width=2.55,height=2.8,pointsize=1)
plot(win,xlim=c(127,143),ylim=c(31,43),col=NA,lwd=0.5)
rect(xleft=100,xright=150,ybottom=20,ytop=50,border=NA,col='powderblue')
abline(v=seq(129,143,2),h=seq(30,42,2),col='grey77',lwd=0.2)
plot(win,xlim=c(127,143),ylim=c(31,43),col='white',add=TRUE,lwd=0.5)
plot(japan,col=adjustcolor(cols,alpha.f = 0.9),border='grey77',add=TRUE,lwd=0.2) 
plot(win,xlim=c(127,143),ylim=c(31,43),col=NA,add=TRUE,lwd=0.5)
plot(c14sites,pch=20,cex=0.2,add=TRUE)
text(x=c(129.99,132.12,131.83,134.91,136.32,141.76,142.58,142.25),y=c(34.24,31.75,35.95,36.42,37.94,35.42,38.54,41.34),labels=as.roman(1:8),cex=1)
axis(1,at=seq(129,143,2),tck=-0.01,padj=-0.8)
axis(2,at=seq(30,42,2),tck=-0.01,padj=0.8)
mtext(side=1,line =1.4,'Longitude')
mtext(side=2,line=1.4,'Latitude')
box()
dev.off()

# Figure 2 ----
load(here('data','sitedata.RData'))
load(here('data','c14data.RData'))

# From Crema et al 2022 (10.1126/sciadv.adc9171), table S4, model b
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
# Figure 3 and 4 ----
load(here('results','icar_c14doubleRes500.RData'))
load(here('results','icar_c14doubleRes750.RData'))
load(here('data','c14data.RData'))


cal500  <- calibrate(c14db500$C14Age,c14db500$C14Error)
cal750  <- calibrate(c14db750$C14Age,c14db750$C14Error)

spds500  <- spds750  <- vector('list',length=8)
riceregions  <- as.character(as.roman(1:8))

for (i in 1:8)
{
	ii  <- which(c14db500$RiceRegion==riceregions[i])
	spds500[[i]]  <- spd(cal500[ii],timeRange=c(c14db500[ii[1],'a'],c14db500[ii[1],'b']),spdnormalised=T)
	ii  <- which(c14db750$RiceRegion==riceregions[i])
	spds750[[i]]  <- spd(cal750[ii],timeRange=c(c14db750[ii[1],'a'],c14db750[ii[1],'b']),spdnormalised=T)
}

pdf(file=here('figures','figure3.pdf'),width=5.3,height=2.5,pointsize=2)
options(scipen=9999)
cex=1
par(mfrow=c(2,4),mar=c(5,3,3,0.5))
for (i in 1:8)
{
	a = c14db500$a[which(c14db500$RiceRegion==as.character(as.roman(i)))][1]
	b = c14db500$b[which(c14db500$RiceRegion==as.character(as.roman(i)))][1]
	deltaT0 = c14db500$ricearrival[which(c14db500$RiceRegion==as.character(as.roman(i)))][1]
	spdBP = spds500[[i]][[2]][,1]
	spdDens1 = spds500[[i]][[2]][,2]

	plot(NULL,xlim=rev(range(spdBP)),ylim=c(0,max(c(spdDens1))),axes=FALSE,xlab='BC/AD',ylab='',cex=cex,cex.lab=cex,cex.main=cex,cex.axis=cex)
	polygon(c(spdBP,rev(spdBP)),c(rep(0,length(spdDens1)),rev(spdDens1)),border=NA,col='lightgrey',cex=cex)
# 	lines(spdBP,spdDens1,lty=1,lwd=0.5)
	abline(v=deltaT0,lty=4,lwd=1)
	modelPlot(model=dDoubleExponentialGrowth,a=a,b=b,params=list(r1=icar.c14double500[,paste0('s1[',i,']')],r2=icar.c14double500[,paste0('s2[',i,']')],mu=rep(deltaT0,nrow(icar.c14double500))),nsample = 100,type='envelope',add=TRUE,col='lightblue',alpha=0.7,lwd=0.7)
	prettyBCAD  <- pretty(BPtoBCAD(spdBP))
	prettyBCAD[which(prettyBCAD==0)] <- 1
	axis(1,at=BCADtoBP(prettyBCAD),labels=abs(prettyBCAD),cex.lab=cex*0.9,cex.axis=cex*0.9,cex=cex*0.85)
	axis(1,at=BCADtoBP(c(seq(-2000,-100,100),1,seq(100,1500,100))),labels=NA,tck=-0.01,cex.axis=cex,cex=cex,cex.lab=cex)
	axis(2,cex.lab=cex,cex.axis=cex*0.9)
	axis(3,at=c(deltaT0-500,deltaT0-250,deltaT0,deltaT0+250,deltaT0+500),labels=c(500,250,0,-250,-500),cex.lab=cex*0.6,cex=cex*0.6,cex.axis=cex*0.9,padj=0.5,gap.axis=-999)
	mtext(TeX(r'($\Delta years$)'),side=3,line=2,cex=cex*0.7,cex.lab=cex*0.7)
	mtext('Summed Probability',side=2,line=2,cex.lab=cex*0.7,cex=cex*0.7)
	legend('topleft',legend=paste0(as.character(as.roman(i))),bty='n',cex=cex*1.5)
	box()
}
dev.off()



pdf(file=here('figures','figure4.pdf'),width=5.3,height=2.5,pointsize=2)
options(scipen=9999)
cex=1
par(mfrow=c(2,4),mar=c(5,3,3,0.5))
for (i in 1:8)
{
	a = c14db750$a[which(c14db750$RiceRegion==as.character(as.roman(i)))][1]
	b = c14db750$b[which(c14db750$RiceRegion==as.character(as.roman(i)))][1]
	deltaT0 = c14db750$ricearrival[which(c14db750$RiceRegion==as.character(as.roman(i)))][1]
	spdBP = spds750[[i]][[2]][,1]
	spdDens1 = spds750[[i]][[2]][,2]

	plot(NULL,xlim=rev(range(spdBP)),ylim=c(0,max(c(spdDens1))),axes=FALSE,xlab='BC/AD',ylab='',cex=cex,cex.lab=cex,cex.main=cex,cex.axis=cex)
	polygon(c(spdBP,rev(spdBP)),c(rep(0,length(spdDens1)),rev(spdDens1)),border=NA,col='lightgrey',cex=cex)
# 	lines(spdBP,spdDens1,lty=1,lwd=0.5)
	abline(v=deltaT0,lty=4,lwd=1)
	modelPlot(model=dDoubleExponentialGrowth,a=a,b=b,params=list(r1=icar.c14double750[,paste0('s1[',i,']')],r2=icar.c14double750[,paste0('s2[',i,']')],mu=rep(deltaT0,nrow(icar.c14double750))),nsample = 100,type='envelope',add=TRUE,col='lightblue',alpha=0.7,lwd=0.7)
	prettyBCAD  <- pretty(BPtoBCAD(spdBP))
	prettyBCAD[which(prettyBCAD==0)] <- 1
	axis(1,at=BCADtoBP(prettyBCAD),labels=abs(prettyBCAD),cex.lab=cex*0.9,cex.axis=cex*0.9,cex=cex*0.85)
	axis(1,at=BCADtoBP(c(seq(-2000,-100,100),1,seq(100,1500,100))),labels=NA,tck=-0.01,cex.axis=cex,cex=cex,cex.lab=cex)
	axis(2,cex.lab=cex,cex.axis=cex*0.9)
	axis(3,at=c(deltaT0-750,deltaT0-500,deltaT0-250,deltaT0,deltaT0+250,deltaT0+500,deltaT0+750),labels=c(750,500,250,0,-250,-500,-750),cex.lab=cex*0.6,cex=cex*0.6,cex.axis=cex*0.9,padj=0.5,gap.axis=-999)
	mtext(TeX(r'($\Delta years$)'),side=3,line=2,cex=cex*0.7,cex.lab=cex*0.7)
	mtext('Summed Probability',side=2,line=2,cex.lab=cex*0.7,cex=cex*0.7)
	legend('topleft',legend=paste0(as.character(as.roman(i))),bty='n',cex=cex*1.5)
	box()
}
dev.off()

# Figure 5 ----
load(here('data','c14data.RData'))
load(here('results','icar_c14doubleRes500.RData'))
load(here('results','icar_c14doubleRes750.RData'))


icar.res.c14.500  <-  icar.c14double500[,1:16] * 100
icar.res.c14.750  <-  icar.c14double750[,1:16] * 100
res500  <- apply(icar.res.c14.500,2,quantile,c(0.05,0.95))
res750  <- apply(icar.res.c14.750,2,quantile,c(0.05,0.95))

ricearrival = data.frame(regions=1:8,m = BCADtoBP(c(-1039,-570,-910,-824,-648,-271,-152,-428)))


main.col = brewer.pal('Set2',n=8)
yrng = range(c(res500,res750))



pdf(file=here('figures','figure5.pdf'),width=5.3,height=5.3,pointsize=7)
par(mfrow=c(2,2),mar=c(2,4,3,1))
plot(NULL,xlim=c(max(ricearrival$m)+800,min(ricearrival$m)-800),ylim=yrng,xlab='',ylab='Annual Growth Rate (%)',axes=F,main='Growth rate 500yrs before farming')
abline(h=seq(-0.2,0.5,0.1),col='lightgrey',lwd=.1)
abline(v=BCADtoBP(c(-1750,-1250,-750,-250,250)),col='lightgrey',lwd=.1)
axis(2,las=2)
box()
axis(1,at=BCADtoBP(c(-1750,-1250,-750,-250,250)),labels=c('1750BC','1250BC','750BC','250BC','250AD'))
abline(h=0,lty=2)

for(i in 1:8)
{
	rect(xleft=ricearrival$m[i]+500,xright=ricearrival$m[i],ybottom=res500[1,i],ytop=res500[2,i],border=main.col[i],lwd=1.5,col=adjustcolor(main.col[i],alpha.f = 0.4))
}


legend('topright',legend=as.roman(1:8),title='Area',fill=adjustcolor(main.col,alpha.f = 0.5),cex=1.4,box.lwd = 1,box.col = 1,bg='white')

plot(NULL,xlim=c(max(ricearrival$m)+800,min(ricearrival$m)-800),ylim=yrng,xlab='',ylab='Annual Growth Rate (%)',axes=F,main='Growth rate 500yrs after farming')
abline(h=seq(-0.2,0.5,0.1),col='lightgrey',lwd=.1)
abline(v=BCADtoBP(c(-1750,-1250,-750,-250,250)),col='lightgrey',lwd=.1)
axis(2,las=2)
box()
axis(1,at=BCADtoBP(c(-1750,-1250,-750,-250,250)),labels=c('1750BC','1250BC','750BC','250BC','250AD'))
abline(h=0,lty=2)
for(i in 1:8)
{
	rect(xleft=ricearrival$m[i],xright=ricearrival$m[i]-500,ybottom=res500[1,i+8],ytop=res500[2,i+8],border=main.col[i],lwd=1.5,col=adjustcolor(main.col[i],alpha.f = 0.4))
}


plot(NULL,xlim=c(max(ricearrival$m)+800,min(ricearrival$m)-800),ylim=yrng,xlab='',ylab='Annual Growth Rate (%)',axes=F,main='Growth rate 750yrs before farming')
abline(h=seq(-0.2,0.5,0.1),col='lightgrey',lwd=.1)
abline(v=BCADtoBP(c(-1750,-1250,-750,-250,250)),col='lightgrey',lwd=.1)
axis(2,las=2)
box()
axis(1,at=BCADtoBP(c(-1750,-1250,-750,-250,250)),labels=c('1750BC','1250BC','750BC','250BC','250AD'))
abline(h=0,lty=2)

for(i in 1:8)
{
	rect(xleft=ricearrival$m[i]+750,xright=ricearrival$m[i],ybottom=res750[1,i],ytop=res750[2,i],border=main.col[i],lwd=1.5,col=adjustcolor(main.col[i],alpha.f = 0.4))
}

plot(NULL,xlim=c(max(ricearrival$m)+800,min(ricearrival$m)-800),ylim=yrng,xlab='',ylab='Annual Growth Rate (%)',axes=F,main='Growth rate 750yrs after farming')

abline(h=seq(-0.2,0.5,0.1),col='lightgrey',lwd=.1)
abline(v=BCADtoBP(c(-1750,-1250,-750,-250,250)),col='lightgrey',lwd=.1)
axis(2,las=2)
box()
axis(1,at=BCADtoBP(c(-1750,-1250,-750,-250,250)),labels=c('1750BC','1250BC','750BC','250BC','250AD'))
abline(h=0,lty=2)

for(i in 1:8)
{
	rect(xleft=ricearrival$m[i],xright=ricearrival$m[i]-750,ybottom=res750[1,i+8],ytop=res750[2,i+8],border=main.col[i],lwd=1.5,col=adjustcolor(main.col[i],alpha.f = 0.4))
}
dev.off()

# Figure 6 ----
load(here('data','c14data.RData'))
load(here('results','icar_c14doubleRes500.RData'))
load(here('results','icar_c14doubleRes750.RData'))

post.bar <- function(x,i,h=0.4,col)
{
	require(grDevices)
	blocks  <- quantile(x,prob=c(0.025,0.1,0.25,0.75,0.9,0.975))
	rect(xleft=blocks[1],xright=blocks[2],ybottom=i-h/5,ytop=i+h/5,border=NA,col=adjustcolor(col,alpha.f=0.3))
	rect(xleft=blocks[2],xright=blocks[3],ybottom=i-h/5,ytop=i+h/5,border=NA,col=adjustcolor(col,alpha.f=0.6))
	rect(xleft=blocks[3],xright=blocks[4],ybottom=i-h/5,ytop=i+h/5,border=NA,col=adjustcolor(col,alpha.f=0.9))
	rect(xleft=blocks[4],xright=blocks[5],ybottom=i-h/5,ytop=i+h/5,border=NA,col=adjustcolor(col,alpha.f=0.6))
	rect(xleft=blocks[5],xright=blocks[6],ybottom=i-h/5,ytop=i+h/5,border=NA,col=adjustcolor(col,alpha.f=0.3))
}




icar.res.c14.500  <- icar.c14double500[,1:16] * 100
icar.res.c14.750  <- icar.c14double750[,1:16] * 100

delta500 <- delta750 <- vector('list',8)

for (i in 1:8)
{
	delta500[[i]] <- icar.res.c14.500[,i+8] - icar.res.c14.500[,i]
	delta750[[i]] <- icar.res.c14.750[,i+8] - icar.res.c14.750[,i]
}




pdf(here('figures','figure6.pdf'),width=2.55,height=3,pointsize=2)
col1  <- "#F95700FF"
col2  <- "#00A4CCFF"
plot(NULL,xlim=c(-0.4,0.8),ylim=c(0.5,23.5),xlab='',ylab='',axes=F)
iseq.a = seq(2,by=3,length.out=8)
iseq.b = seq(1,by=3,length.out=8)
abline(h=seq(3,by=3,length.out=7),col='darkgrey',lty=2)

for (i in 1:8)
{
	b  <- delta500[[i]]
	a  <- delta750[[i]]
	post.bar(b,i=iseq.b[i],col=col1,h=1.5)
	post.bar(a,i=iseq.a[i],col=col2,h=1.5)
}

axis(2,at=iseq.a-0.5,labels = c('I','II','III','IV','V','VI','VII','VIII'),las=2)
axis(1)
abline(v=0,lty=3)
rect(xleft=0.43,xright=0.85,ybottom=14,ytop=24.5,border=NA,col='white')
legend(x=0.45,y=24,legend=c('95% HPD','80% HPD','50% HPD'),fill=c(adjustcolor(col2,alpha.f = 0.9),adjustcolor(col2,alpha.f = 0.6),adjustcolor(col2,alpha.f = 0.3)),title='750 yrs window',bty='n',bg='white',cex=0.9)
legend(x=0.45,y=19,legend=c('95% HPD','80% HPD','50% HPD'),fill=c(adjustcolor(col1,alpha.f = 0.9),adjustcolor(col1,alpha.f = 0.6),adjustcolor(col1,alpha.f = 0.3)),title='500 yrs window',bty='n',bg='white',cex=0.9)
box()
mtext(text=TeX(r'($\Delta r)'),side=1,line=2.5)
mtext(text='Area',side=2,line=2.5)
dev.off()
