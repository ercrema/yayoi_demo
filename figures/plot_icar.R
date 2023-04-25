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

# Load ICAR Results ----
load(here('results','icar_c14doubleRes.RData'))
load(here('results','icar_c14doubleRes750.RData'))
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



# Summarise data ----
icar.res.c14.500  <- icar.res.c14.500.net <-  icar.c14double[,1:16] * 100
icar.res.c14.750  <- icar.res.c14.750.net <-  icar.c14double750[,1:16] * 100
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


# Plot Comparison ----
main.col <- c(rgb(51,34,136,maxColorValue=255),rgb(136,204,238,maxColorValue = 255),rgb(68,170,153,maxColorValue = 255),rgb(17,119,51,maxColorValue = 255),rgb(153,153,51,maxColorValue = 255),rgb(221,204,119,maxColorValue = 255),rgb(204,102,119,maxColorValue = 255),rgb(136,34,85,maxColorValue = 255))


win  <- ne_countries(continent = 'asia',scale=10)
japan <- ne_states(country = "japan") |> subset(!name_vi %in% c("Okinawa","Hokkaido"))
df.pref.reg <- read.csv(here('data','prefecture_region_match.csv'))
japan@data <- left_join(japan@data,df.pref.reg,by=c('name'='Prefecture'))
japan <- gUnaryUnion(japan,id=japan@data$Area)
japan.sf <- as(japan,'sf')
win <- gUnaryUnion(win,id=win@data[,1])


par(mfrow=c(1,3))
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
