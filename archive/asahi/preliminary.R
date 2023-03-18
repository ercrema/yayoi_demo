# Load Libray & Data ----
library(here)
library(rcarbon)
load(here('data','c14data.RData'))
load(here('data','rekihaku_jomonyayoi_settlement_db_v01.RData'))

# Subset to radius around Asahi-site ----
asahi  <- c(35.221111,136.854444)
asahi.dates  <- subset(c14db, SiteName=='朝日遺跡' & PrefectureCode==55 & Material!='Marine')
cal.asahi.dates <- calibrate(asahi.dates$C14Age,asahi.dates$C14Error)
asahi.all  <- spd(cal.asahi.dates,timeRange=c(3500,1500),runm=10)
asahi.anthropic.nonora  <- spd(cal.asahi.dates[which(asahi.dates$anthropic==TRUE & asahi.dates$Material=='Terrestrial')],timeRange=c(3500,1500),runm=10)
asahi.ora  <- spd(cal.asahi.dates[which(asahi.dates$Material_Details=='Organic Residue from Pottery')],timeRange=c(3500,1500),runm=10)
pdf(file=here('asahi','asahi_dates.pdf'),width=6,height=8)
par(mfrow=c(3,1))
plot(asahi.all,calendar='BCAD',main='Asahi Site - All Dates',ylim=c(0,0.18))
arrows(x0=-771,x1=-552,y0=0.1,y1=0.1,length = 0.05,code=3,angle=90)
text(x=median(c(-771,-552)),y=0.135,label='Estimated Arrival of \n Rice Farming (Region V)')
legend('topleft',legend=c(paste('n=',nrow(asahi.dates))),bty='n',cex=2)
barCodes(x=BPtoBCAD(medCal(cal.asahi.dates)),y=rev(quantile(c(par('usr')[3:4]),prob=c(0.95,1))),width = 5,col = rgb(0,0,0,0.3))
plot(asahi.anthropic.nonora,calendar='BCAD',fill='lightblue',main = 'Asahi Site - Anthropogeneic Dates (excluding ORA)',ylim=c(0,0.18))
legend('topleft',legend=c(paste('n=',sum(asahi.dates$anthropic==TRUE & asahi.dates$Material=='Terrestrial'))),bty='n',cex=2)
barCodes(x=BPtoBCAD(medCal(cal.asahi.dates[which(asahi.dates$anthropic==TRUE & asahi.dates$Material=='Terrestrial')])),y=rev(quantile(c(par('usr')[3:4]),prob=c(0.95,1))),width = 5,col = rgb(0,0,0,0.3))
plot(asahi.ora,calendar='BCAD',fill='darkred',main = 'Asahi Site - Dates from OR',ylim=c(0,0.18))
legend('topleft',legend=c(paste('n=',sum(asahi.dates$Material_Details=='Organic Residue from Pottery'))),bty='n',cex=2)
barCodes(x=BPtoBCAD(medCal(cal.asahi.dates[which(asahi.dates$Material_Details=='Organic Residue from Pottery')])),y=rev(quantile(c(par('usr')[3:4]),prob=c(0.95,1))),width = 5,col = rgb(0,0,0,0.3))
dev.off()


