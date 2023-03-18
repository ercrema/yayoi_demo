library(here)
library(nimbleCarbon)
library(RColorBrewer)
library(dplyr)
library(rcarbon)

# Define Region Colours ----
main.col <- brewer.pal('Set3',n=7)

# Obtain Sample Sizes ----
# Load Data ----
load(here('data','c14data.RData'))
c14db = subset(c14db,C14Age<4000 & C14Age> 300 &  Material == 'Terrestrial' & anthropic == TRUE) |> select(C14Age,C14Error,SiteID,Region=Region3,RegionID=RegionID3) |> arrange(RegionID)


# Site Level Thinning ----
calibrated.dates <- calibrate(c14db$C14Age,c14db$C14Error)
bin  <- binPrep(sites = c14db$SiteID,ages=calibrated.dates,h=50)
i  <- thinDates(ages=c14db$C14Age,  errors=c14db$C14Error, bins=bin, size=1, thresh=1,seed=123,method='splitsample')
c14db  <-  c14db[i,]
calibrated.dates <- calibrated.dates[i]

# Setup Regions ----
gapYear <- 800 
regionStartDates = unique(select(c14db,Region,RegionID)) |> arrange(RegionID)
regionStartDates$a1 = BCADtoBP(c(-1283,-1117,-981,-833,-685,-515,-362)) # 1283, 1117, 981, 833, 685, 515, 362
regionStartDates$a2 = BCADtoBP(c(-900))
regionStartDates$b1 = regionStartDates$a1 - gapYear
regionStartDates$b2 = BCADtoBP(250)

# Subset Data to Relevant Temporal Window of Analyses ----
c14db$med = medCal(calibrated.dates)
c14db$r1 = NA
c14db$r2 = NA

for (i in 1:nrow(regionStartDates))
{
	a1 = regionStartDates$a1[i]
	a2 = regionStartDates$a2[i]
	b1 = regionStartDates$b1[i]
	b2 = regionStartDates$b2[i]
	ii = which(c14db$RegionID == i)
	c14db$r1[ii] = FALSE
	c14db$r2[ii] = FALSE
	regionalDates = calibrated.dates[ii]
	iii = ii[which.CalDates(regionalDates,BP < a1 & BP > b1,p=0.5)]
	c14db$r1[iii] = TRUE
	iii = ii[which.CalDates(regionalDates,BP < a2 & BP > b2,p=0.5)]
	c14db$r2[iii] = TRUE
}

# Define Data & Constants ----
data(intcal20)
c14db1 = subset(c14db,r1==TRUE)
c14db2 = subset(c14db,r2==TRUE)

nsamples_800yrs <- table(c14db1$RegionID)
nsamples_yayoi <- table(c14db2$RegionID)



# Load Results
load(here('R_images','regdemo_800yrs.RData'))
load(here('R_images','regdemo_yayoi.RData'))

Regions = c('Kyushu','Chugoku + Shikoku','Kansai','Chubu','Kanto','S. Tohoku + Niigata','N. Tohoku')

pdf(file = here('figures','regdemoyayoi_results.pdf'),width = 7, height=12)
index <- c(7,6,5,4,3,1) #Excluding Chugoku/Shikoku

par(mfrow=c(6,1),mar=c(4,1,0,0))
for (i in index)
{
  samples = post.sample.yayoi[[1]][,i] * 100
  ndates = nsamples_yayoi[[i]]
  plottitle = paste0(Regions[i],'\nNumber of Dates = ',ndates)
  postHPDplot(samples,xlim=c(-0.2,0.35),xlab='Annual Growth Rate (%)',main='',ylab='Probability Density',axes=F,hpd.col=main.col[i])
  abline(v=0,lty=2)
  axis(1)
  legend('topleft',legend=plottitle,bty='n',cex=1.4)
}
dev.off()



pdf(file = here('figures','regdemo800yrs_results.pdf'),width = 7, height=12)
index <- c(7,6,5,4,3,1) #Excluding Chugoku/Shikoku

par(mfrow=c(6,1),mar=c(4,1,0,0))
for (i in index)
{
  samples = post.sample.800yrs[[1]][,i] * 100
  ndates = nsamples_800yrs[[i]]
  plottitle = paste0(Regions[i],'\nNumber of Dates = ',ndates)
  postHPDplot(samples,xlim=c(-0.2,0.35),xlab='Annual Growth Rate (%)',main='',ylab='Probability Density',axes=F,hpd.col=main.col[i])
  abline(v=0,lty=2)
  axis(1)
  legend('topleft',legend=plottitle,bty='n',cex=1.4)
}
dev.off()

# Posterior Predictive Check ----
load(here('R_images','ppcheck_800yrs.RData'))
load(here('R_images','ppcheck_yayoi.RData'))


pdf(file = here('figures','ppcheck800yrs_results.pdf'),width = 6, height=12)
par(mfrow=c(6,1),mar=c(4.3,3,2,1))
for (i in index)
{
  ndates = nsamples_800yrs[[i]]
  cc <- postPredCor(ppcheck_800yrs[[i]])
  cc.med <- round(median(cc),2)
  cc.lo <- round(quantile(cc,0.025),2)
  cc.hi <- round(quantile(cc,0.975),2)
  plottitle = paste0(Regions[i],'\nNumber of Dates: ',ndates,'\nGoodness of fit: ',cc.med,' (',cc.lo,'~',cc.hi,')')
  plot(ppcheck_800yrs[[i]],calendar='BCAD')
  yrs = c(1,200,400,600,800,1000)
  ats = BPtoBCAD((regionStartDates$a1[i]:regionStartDates$b1[i])[yrs])
  axis(3,at = ats,labels = yrs,padj = 1.2,cex=0.8,tck=-.01)
#   mtext('Years from Rice Farming',side=3,line=1.5,cex=0.5)
  legend('topleft',legend=plottitle,bty='n',cex=1)
}
dev.off() 


