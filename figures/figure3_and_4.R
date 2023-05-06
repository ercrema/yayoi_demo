library(here)
library(latex2exp)
library(nimbleCarbon)
library(rcarbon)

# Load ICAR Results ----
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


