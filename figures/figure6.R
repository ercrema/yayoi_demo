# Load Library and Data ----
library(here)
library(rcarbon)
library(latex2exp)
library(RColorBrewer)
library(grDevices)
load(here('data','c14data.RData'))
load(here('results','icar_c14doubleRes500.RData'))
load(here('results','icar_c14doubleRes750.RData'))


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
legend('topright',legend=c('750yrs window','500yrs window'),fill=c(col2,col1),bty='n')
box()
mtext(text=TeX(r'($\Delta r)'),side=1,line=2.5)
mtext(text='Area',side=2,line=2.5)
dev.off()
