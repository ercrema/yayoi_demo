library(ggplot2)
library(dplyr)
library(gridExtra)
library(tidyr)
library(grid)
library(nimbleCarbon)
library(rcarbon)
library(rnaturalearth)
library(sf)
library(maptools)
library(rgeos)
library(tidybayes)
library(spdep)
library(here)
library(parallel)
library(latex2exp)

# Load ICAR Results ----
load(here('results','icar_c14doubleRes500.RData'))
load(here('results','icar_c14doubleRes750.RData'))
load(here('results','icar_c14doubleRes1000.RData'))
load(here('results','icar_sitesdoubleRes.RData'))


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

cal500  <- calibrate(c14db500$C14Age,c14db500$C14Error)
cal750  <- calibrate(c14db750$C14Age,c14db750$C14Error)
cal1000  <- calibrate(c14db1000$C14Age,c14db1000$C14Error)

cal500n  <- calibrate(c14db500$C14Age,c14db500$C14Error,normalised=F)
cal750n  <- calibrate(c14db750$C14Age,c14db750$C14Error,normalised=F)
cal1000n  <- calibrate(c14db1000$C14Age,c14db1000$C14Error,normalised=F)

spds500  <- spds750  <- spds1000 <- vector('list',length=8)
spds500n  <- spds750n  <- spds1000n <- vector('list',length=8)
riceregions  <- as.character(as.roman(1:8))

for (i in 1:8)
{
	ii  <- which(c14db500$RiceRegion==riceregions[i])
	spds500[[i]]  <- spd(cal500[ii],timeRange=c(c14db500[ii[1],'a'],c14db500[ii[1],'b']),spdnormalised=T)
	spds500n[[i]]  <- spd(cal500n[ii],timeRange=c(c14db500[ii[1],'a'],c14db500[ii[1],'b']),spdnormalised=T)

	ii  <- which(c14db750$RiceRegion==riceregions[i])
	spds750[[i]]  <- spd(cal750[ii],timeRange=c(c14db750[ii[1],'a'],c14db750[ii[1],'b']),spdnormalised=T)
	spds750n[[i]]  <- spd(cal750n[ii],timeRange=c(c14db750[ii[1],'a'],c14db750[ii[1],'b']),spdnormalised=T)

	ii  <- which(c14db1000$RiceRegion==riceregions[i])
	spds1000[[i]]  <- spd(cal1000[ii],timeRange=c(c14db1000[ii[1],'a'],c14db1000[ii[1],'b']),spdnormalised=T)
	spds1000n[[i]]  <- spd(cal1000n[ii],timeRange=c(c14db1000[ii[1],'a'],c14db1000[ii[1],'b']),spdnormalised=T)
}

# par(mfrow=c(2,4))
# lapply(spds500n,plot)
# lapply(spds750n,plot)
# lapply(spds1000n,plot,calendar='BCAD')

par(mfrow=c(2,4))
for (i in 1:8)
{
	a = c14db1000$a[which(c14db1000$RiceRegion==as.character(as.roman(i)))][1]
	b = c14db1000$b[which(c14db1000$RiceRegion==as.character(as.roman(i)))][1]
	deltaT0 = c14db1000$ricearrival[which(c14db1000$RiceRegion==as.character(as.roman(i)))][1]
	spdBP = spds1000[[i]][[2]][,1]
	spdDens1 = spds1000[[i]][[2]][,2]
	spdDens2 = spds1000n[[i]][[2]][,2]

	plot(NULL,xlim=rev(range(spdBP)),ylim=c(0,max(c(spdDens1,spdDens2))),axes=FALSE,xlab='BC/AD',ylab='Summed Probability')
	polygon(c(spdBP,rev(spdBP)),c(rep(0,length(spdDens2)),rev(spdDens2)),border=NA,col='lightgrey')
	lines(spdBP,spdDens1,lty=2)
	abline(v=deltaT0,lty=4,lwd=1)
	modelPlot(model=dDoubleExponentialGrowth,a=a,b=b,params=list(r1=icar.c14double1000[,paste0('s1[',i,']')],r2=icar.c14double1000[,paste0('s2[',i,']')],mu=rep(deltaT0,nrow(icar.c14double1000))),nsample = 100,type='envelope',add=TRUE,col='lightblue',alpha=0.5)
	prettyBCAD  <- pretty(BPtoBCAD(spdBP))
	prettyBCAD[which(prettyBCAD==0)] <- 1
	axis(1,at=BCADtoBP(prettyBCAD),labels=abs(prettyBCAD))
	axis(1,at=BCADtoBP(c(seq(-2000,-100,100),1,seq(100,1500,100))),labels=NA,tck=-0.01)
	axis(2)
	axis(3,at=c(deltaT0-1000,deltaT0-750,deltaT0-500,deltaT0-250,deltaT0,deltaT0+250,deltaT0+500,deltaT0+750,deltaT0+1000),labels=c(-1000,-750,-500,-250,0,250,500,750,1000))
	mtext(TeX(r'($\Delta years$)'),side=3,line=2)
	legend('topleft',legend=paste0('Area ',as.character(as.roman(i))),cex=1.5,bty='n')
	box()
}

par(mfrow=c(2,4))
for (i in 1:8)
{
	a = c14db750$a[which(c14db750$RiceRegion==as.character(as.roman(i)))][1]
	b = c14db750$b[which(c14db750$RiceRegion==as.character(as.roman(i)))][1]
	deltaT0 = c14db750$ricearrival[which(c14db750$RiceRegion==as.character(as.roman(i)))][1]
	spdBP = spds750[[i]][[2]][,1]
	spdDens1 = spds750[[i]][[2]][,2]
	spdDens2 = spds750n[[i]][[2]][,2]

	plot(NULL,xlim=rev(range(spdBP)),ylim=c(0,max(c(spdDens1,spdDens2))),axes=FALSE,xlab='BC/AD',ylab='Summed Probability')
	polygon(c(spdBP,rev(spdBP)),c(rep(0,length(spdDens2)),rev(spdDens2)),border=NA,col='lightgrey')
	lines(spdBP,spdDens1,lty=2)
	abline(v=deltaT0,lty=4,lwd=1)
	modelPlot(model=dDoubleExponentialGrowth,a=a,b=b,params=list(r1=icar.c14double750[,paste0('s1[',i,']')],r2=icar.c14double750[,paste0('s2[',i,']')],mu=rep(deltaT0,nrow(icar.c14double750))),nsample = 100,type='envelope',add=TRUE,col='lightblue',alpha=0.5)
	prettyBCAD  <- pretty(BPtoBCAD(spdBP))
	prettyBCAD[which(prettyBCAD==0)] <- 1
	axis(1,at=BCADtoBP(prettyBCAD),labels=abs(prettyBCAD))
	axis(1,at=BCADtoBP(c(seq(-2000,-100,100),1,seq(100,1500,100))),labels=NA,tck=-0.01)
	axis(2)
	axis(3,at=c(deltaT0-750,deltaT0-500,deltaT0-250,deltaT0,deltaT0+250,deltaT0+500,deltaT0+750),labels=c(-750,-500,-250,0,250,500,750))
	mtext(TeX(r'($\Delta years$)'),side=3,line=2)
	legend('topleft',legend=paste0('Area ',as.character(as.roman(i))),cex=1.5,bty='n')
	box()
}



par(mfrow=c(2,4))
for (i in 1:8)
{
	a = c14db500$a[which(c14db500$RiceRegion==as.character(as.roman(i)))][1]
	b = c14db500$b[which(c14db500$RiceRegion==as.character(as.roman(i)))][1]
	deltaT0 = c14db500$ricearrival[which(c14db500$RiceRegion==as.character(as.roman(i)))][1]
	spdBP = spds500[[i]][[2]][,1]
	spdDens1 = spds500[[i]][[2]][,2]
	spdDens2 = spds500n[[i]][[2]][,2]

	plot(NULL,xlim=rev(range(spdBP)),ylim=c(0,max(c(spdDens1,spdDens2))),axes=FALSE,xlab='BC/AD',ylab='Summed Probability')
	polygon(c(spdBP,rev(spdBP)),c(rep(0,length(spdDens2)),rev(spdDens2)),border=NA,col='lightgrey')
	lines(spdBP,spdDens1,lty=2)
	abline(v=deltaT0,lty=4,lwd=1)
	modelPlot(model=dDoubleExponentialGrowth,a=a,b=b,params=list(r1=icar.c14double500[,paste0('s1[',i,']')],r2=icar.c14double500[,paste0('s2[',i,']')],mu=rep(deltaT0,nrow(icar.c14double500))),nsample = 100,type='envelope',add=TRUE,col='lightblue',alpha=0.5)
	prettyBCAD  <- pretty(BPtoBCAD(spdBP))
	prettyBCAD[which(prettyBCAD==0)] <- 1
	axis(1,at=BCADtoBP(prettyBCAD),labels=abs(prettyBCAD))
	axis(1,at=BCADtoBP(c(seq(-2000,-100,100),1,seq(100,1500,100))),labels=NA,tck=-0.01)
	axis(2)
	axis(3,at=c(deltaT0-500,deltaT0-250,deltaT0,deltaT0+250,deltaT0+500),labels=c(-500,-250,0,250,500))
	mtext(TeX(r'($\Delta years$)'),side=3,line=2)
	legend('topleft',legend=paste0('Area ',as.character(as.roman(i))),cex=1.5,bty='n')
	box()
}


# Summarise data ----
icar.res.c14.500  <- icar.res.c14.500.net <-  icar.c14double[,1:16] * 100
icar.res.c14.750  <- icar.res.c14.750.net <-  icar.c14double750[,1:16] * 100
icar.res.c14.1000  <- icar.res.c14.1000.net <-  icar.c14double1000[,1:16] * 100
icar.res.sites.500  <- icar.res.sites.500.net <-  icar.sitesdouble[,1:16] * 100

icar.res.c14.500.net$delta1  <- icar.res.c14.500[,'s2[1]'] - icar.res.c14.500[,'s1[1]']
icar.res.c14.500.net$delta2  <- icar.res.c14.500[,'s2[2]'] - icar.res.c14.500[,'s1[2]']
icar.res.c14.500.net$delta3  <- icar.res.c14.500[,'s2[3]'] - icar.res.c14.500[,'s1[3]']
icar.res.c14.500.net$delta4  <- icar.res.c14.500[,'s2[4]'] - icar.res.c14.500[,'s1[4]']
icar.res.c14.500.net$delta5  <- icar.res.c14.500[,'s2[5]'] - icar.res.c14.500[,'s1[5]']
icar.res.c14.500.net$delta6  <- icar.res.c14.500[,'s2[6]'] - icar.res.c14.500[,'s1[6]']
icar.res.c14.500.net$delta7  <- icar.res.c14.500[,'s2[7]'] - icar.res.c14.500[,'s1[7]']
icar.res.c14.500.net$delta8  <- icar.res.c14.500[,'s2[8]'] - icar.res.c14.500[,'s1[8]']
icar.res.c14.500.net  <- icar.res.c14.500.net[,grep("delta",colnames(icar.res.c14.500.net))]

icar.res.sites.500.net$delta1  <- icar.res.sites.500[,'s2[1]'] - icar.res.sites.500[,'s1[1]']
icar.res.sites.500.net$delta2  <- icar.res.sites.500[,'s2[2]'] - icar.res.sites.500[,'s1[2]']
icar.res.sites.500.net$delta3  <- icar.res.sites.500[,'s2[3]'] - icar.res.sites.500[,'s1[3]']
icar.res.sites.500.net$delta4  <- icar.res.sites.500[,'s2[4]'] - icar.res.sites.500[,'s1[4]']
icar.res.sites.500.net$delta5  <- icar.res.sites.500[,'s2[5]'] - icar.res.sites.500[,'s1[5]']
icar.res.sites.500.net$delta6  <- icar.res.sites.500[,'s2[6]'] - icar.res.sites.500[,'s1[6]']
icar.res.sites.500.net$delta7  <- icar.res.sites.500[,'s2[7]'] - icar.res.sites.500[,'s1[7]']
icar.res.sites.500.net$delta8  <- icar.res.sites.500[,'s2[8]'] - icar.res.sites.500[,'s1[8]']
icar.res.sites.500.net  <- icar.res.sites.500.net[,grep("delta",colnames(icar.res.sites.500.net))]

icar.res.c14.750.net$delta1  <- icar.res.c14.750[,'s2[1]'] - icar.res.c14.750[,'s1[1]']
icar.res.c14.750.net$delta2  <- icar.res.c14.750[,'s2[2]'] - icar.res.c14.750[,'s1[2]']
icar.res.c14.750.net$delta3  <- icar.res.c14.750[,'s2[3]'] - icar.res.c14.750[,'s1[3]']
icar.res.c14.750.net$delta4  <- icar.res.c14.750[,'s2[4]'] - icar.res.c14.750[,'s1[4]']
icar.res.c14.750.net$delta5  <- icar.res.c14.750[,'s2[5]'] - icar.res.c14.750[,'s1[5]']
icar.res.c14.750.net$delta6  <- icar.res.c14.750[,'s2[6]'] - icar.res.c14.750[,'s1[6]']
icar.res.c14.750.net$delta7  <- icar.res.c14.750[,'s2[7]'] - icar.res.c14.750[,'s1[7]']
icar.res.c14.750.net$delta8  <- icar.res.c14.750[,'s2[8]'] - icar.res.c14.750[,'s1[8]']
icar.res.c14.750.net  <- icar.res.c14.750.net[,grep("delta",colnames(icar.res.c14.750.net))]

icar.res.c14.1000.net$delta1  <- icar.res.c14.1000[,'s2[1]'] - icar.res.c14.1000[,'s1[1]']
icar.res.c14.1000.net$delta2  <- icar.res.c14.1000[,'s2[2]'] - icar.res.c14.1000[,'s1[2]']
icar.res.c14.1000.net$delta3  <- icar.res.c14.1000[,'s2[3]'] - icar.res.c14.1000[,'s1[3]']
icar.res.c14.1000.net$delta4  <- icar.res.c14.1000[,'s2[4]'] - icar.res.c14.1000[,'s1[4]']
icar.res.c14.1000.net$delta5  <- icar.res.c14.1000[,'s2[5]'] - icar.res.c14.1000[,'s1[5]']
icar.res.c14.1000.net$delta6  <- icar.res.c14.1000[,'s2[6]'] - icar.res.c14.1000[,'s1[6]']
icar.res.c14.1000.net$delta7  <- icar.res.c14.1000[,'s2[7]'] - icar.res.c14.1000[,'s1[7]']
icar.res.c14.1000.net$delta8  <- icar.res.c14.1000[,'s2[8]'] - icar.res.c14.1000[,'s1[8]']
icar.res.c14.1000.net  <- icar.res.c14.1000.net[,grep("delta",colnames(icar.res.c14.1000.net))]

c14.500  <- gather(icar.res.c14.500)
c14.500$periods  <- NA
c14.500$periods  <- "before"
c14.500$periods[grep("s2",c14.500$key)]  <- "after"
c14.500$region  <- as.roman(as.numeric(substr(c14.500$key,4,4)))

sites.500  <- gather(icar.res.sites.500)
sites.500$periods  <- NA
sites.500$periods  <- "before"
sites.500$periods[grep("s2",sites.500$key)]  <- "after"
sites.500$region  <- as.roman(as.numeric(substr(sites.500$key,4,4)))

c14.750  <- gather(icar.res.c14.750)
c14.750$periods  <- NA
c14.750$periods  <- "before"
c14.750$periods[grep("s2",c14.750$key)]  <- "after"
c14.750$region  <- as.roman(as.numeric(substr(c14.750$key,4,4)))

c14.1000  <- gather(icar.res.c14.1000)
c14.1000$periods  <- NA
c14.1000$periods  <- "before"
c14.1000$periods[grep("s2",c14.1000$key)]  <- "after"
c14.1000$region  <- as.roman(as.numeric(substr(c14.1000$key,4,4)))


# Plot Comparison ----
main.col <- c(rgb(51,34,136,maxColorValue=255),rgb(136,204,238,maxColorValue = 255),rgb(68,170,153,maxColorValue = 255),rgb(17,119,51,maxColorValue = 255),rgb(153,153,51,maxColorValue = 255),rgb(221,204,119,maxColorValue = 255),rgb(204,102,119,maxColorValue = 255),rgb(136,34,85,maxColorValue = 255))


win  <- ne_countries(continent = 'asia',scale=10)
japan <- ne_states(country = "japan") |> subset(!name_vi %in% c("Okinawa","Hokkaido"))
df.pref.reg <- read.csv(here('data','prefecture_region_match.csv'))
japan@data <- left_join(japan@data,df.pref.reg,by=c('name'='Prefecture'))
japan <- gUnaryUnion(japan,id=japan@data$Area)
japan.sf <- as(japan,'sf')
win <- gUnaryUnion(win,id=win@data[,1])


par(mfrow=c(2,2))
plot(win,xlim=c(127,143),ylim=c(31,43),col='lightgrey')
plot(japan.sf,col='darkgrey',border=1,add=TRUE)
text(x=c(129.99,132.12,131.83,134.91,136.32,141.76,142.58,142.25),y=c(34.24,31.75,35.95,36.42,37.94,35.42,38.54,41.34),labels=as.roman(1:8),cex=1.5)
axis(1,at=seq(129,143,2),tck=-0.01,padj=-0.8)
axis(2,at=seq(30,42,2),tck=-0.01,padj=0.8)
mtext(side=1,line =1.4,'Longitude')
mtext(side=2,line=1.4,'Latitude')
box()

plot(NULL,xlim=c(-0.5,0.5),ylim=c(0.5,23.5),xlab='',ylab='',axes=F,main='500yrs before vs after - 14C dates')
iseq.a = seq(2,by=3,length.out=8)
iseq.b = seq(1,by=3,length.out=8)
abline(h=seq(3,by=3,length.out=7),col='darkgrey',lty=2)

for (i in 1:8)
{
	b  <- subset(c14.500,region==i&periods=='before')$value
	a  <- subset(c14.500,region==i&periods=='after')$value
	post.bar(b,i=iseq.b[i],col='royalblue',h=1)
	post.bar(a,i=iseq.a[i],col='indianred',h=1)
}

axis(2,at=iseq.a-0.5,labels = c('I','II','III','IV','V','VI','VII','VIII'),las=2)
axis(1)
abline(v=0,lty=3)
box()

plot(NULL,xlim=c(-0.5,0.5),ylim=c(0.5,23.5),xlab='',ylab='',axes=F,main='1000yrs before vs after - 14C dates')
iseq.a = seq(2,by=3,length.out=8)
iseq.b = seq(1,by=3,length.out=8)
abline(h=seq(3,by=3,length.out=7),col='darkgrey',lty=2)

for (i in 1:8)
{
	b  <- subset(c14.1000,region==i&periods=='before')$value
	a  <- subset(c14.1000,region==i&periods=='after')$value
	post.bar(b,i=iseq.b[i],col='royalblue',h=1)
	post.bar(a,i=iseq.a[i],col='indianred',h=1)
}

axis(2,at=iseq.a-0.5,labels = c('I','II','III','IV','V','VI','VII','VIII'),las=2)
axis(1)
abline(v=0,lty=3)
box()



plot(NULL,xlim=c(-0.5,0.5),ylim=c(0.5,23.5),xlab='',ylab='',axes=F,main='750yrs before vs after - 14C dates')
iseq.a = seq(2,by=3,length.out=8)
iseq.b = seq(1,by=3,length.out=8)
abline(h=seq(3,by=3,length.out=7),col='darkgrey',lty=2)

for (i in 1:8)
{
	b  <- subset(c14.750,region==i&periods=='before')$value
	a  <- subset(c14.750,region==i&periods=='after')$value
	post.bar(b,i=iseq.b[i],col='royalblue',h=1)
	post.bar(a,i=iseq.a[i],col='indianred',h=1)
}

axis(2,at=iseq.a-0.5,labels = c('I','II','III','IV','V','VI','VII','VIII'),las=2)
axis(1)
abline(v=0,lty=3)
box()
