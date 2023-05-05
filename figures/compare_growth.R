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


par(mfrow=c(2,2))
plot(NULL,xlim=c(max(ricearrival$m)+800,min(ricearrival$m)-800),ylim=yrng,xlab='BC/AD',ylab='Annual Growth Rate (%)',axes=F,main='Growth rate 500yrs before farming')
axis(2,las=2)
box()
axis(1,at=BCADtoBP(c(seq(-2000,-250,250),1,c(250,500,750))),labels=c(seq(2000,250,-250),1,250,500,750))
abline(h=0,lty=2)

for(i in 1:8)
{
	rect(xleft=ricearrival$m[i]+500,xright=ricearrival$m[i],ybottom=res500[1,i],ytop=res500[2,i],border=main.col[i],lwd=1.5,col=adjustcolor(main.col[i],alpha.f = 0.4))
}


legend('topright',legend=as.roman(1:8),title='Area',bty='n',fill=adjustcolor(main.col,alpha.f = 0.5),cex=1.2)

plot(NULL,xlim=c(max(ricearrival$m)+800,min(ricearrival$m)-800),ylim=yrng,xlab='BC/AD',ylab='Annual Growth Rate (%)',axes=F,main='Growth rate 500yrs after farming')
axis(2,las=2)
box()
axis(1,at=BCADtoBP(c(seq(-2000,-250,250),1,c(250,500,750))),labels=c(seq(2000,250,-250),1,250,500,750))
abline(h=0,lty=2)

for(i in 1:8)
{
	rect(xleft=ricearrival$m[i],xright=ricearrival$m[i]-500,ybottom=res500[1,i+8],ytop=res500[2,i+8],border=main.col[i],lwd=1.5,col=adjustcolor(main.col[i],alpha.f = 0.4))
}


plot(NULL,xlim=c(max(ricearrival$m)+800,min(ricearrival$m)-800),ylim=yrng,xlab='BC/AD',ylab='Annual Growth Rate (%)',axes=F,main='Growth rate 750yrs before farming')
axis(2,las=2)
box()
axis(1,at=BCADtoBP(c(seq(-2000,-250,250),1,c(250,500,750))),labels=c(seq(2000,250,-250),1,250,500,750))
abline(h=0,lty=2)

for(i in 1:8)
{
	rect(xleft=ricearrival$m[i]+750,xright=ricearrival$m[i],ybottom=res750[1,i],ytop=res750[2,i],border=main.col[i],lwd=1.5,col=adjustcolor(main.col[i],alpha.f = 0.4))
}

plot(NULL,xlim=c(max(ricearrival$m)+800,min(ricearrival$m)-800),ylim=yrng,xlab='BC/AD',ylab='Annual Growth Rate (%)',axes=F,main='Growth rate 750yrs after farming')
axis(2,las=2)
box()
axis(1,at=BCADtoBP(c(seq(-2000,-250,250),1,c(250,500,750))),labels=c(seq(2000,250,-250),1,250,500,750))
abline(h=0,lty=2)

for(i in 1:8)
{
	rect(xleft=ricearrival$m[i],xright=ricearrival$m[i]-750,ybottom=res750[1,i+8],ytop=res750[2,i+8],border=main.col[i],lwd=1.5,col=adjustcolor(main.col[i],alpha.f = 0.4))
}


# Delta Plot ----
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


plot(NULL,xlim=c(-0.4,0.8),ylim=c(0.5,23.5),xlab='',ylab='',axes=F)
iseq.a = seq(2,by=3,length.out=8)
iseq.b = seq(1,by=3,length.out=8)
abline(h=seq(3,by=3,length.out=7),col='darkgrey',lty=2)

for (i in 1:8)
{
	b  <- delta500[[i]]
	a  <- delta750[[i]]
	post.bar(b,i=iseq.b[i],col='darkgreen',h=1.2)
	post.bar(a,i=iseq.a[i],col='darkorange',h=1.2)
}

axis(2,at=iseq.a-0.5,labels = c('I','II','III','IV','V','VI','VII','VIII'),las=2)
axis(1)
abline(v=0,lty=3)
legend('topright',legend=c('500yrs window','750yrs window'),fill=c('darkgreen','darkorange'))
box()
mtext(text=TeX(r'($\Delta r)'),side=1,line=2.5)

