# Load Library and Data ----
library(here)
library(rcarbon)
library(latex2exp)
library(RColorBrewer)
library(grDevices)
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

