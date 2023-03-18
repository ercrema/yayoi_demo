library(here)
library(rcarbon)
# Settings
timeRange = c(4049,1700)


# Load R functions and Data
load(here('data','rekihaku_jomonyayoi_settlement_db_v01.RData'))
source(here('R','ageSubset.R'))
source(here('R','countSites.R'))
source(here('R','countSitesPermTest.R'))
source(here('R','plotPermResult.R'))


# Compute Observed Frequencies ----
## All Regions ----
sitesFreqAllRegions= countSites(ageRange=timeRange,ageMatrix = ageMatrix,sites=sites,ageRef=ages)
pdf(file = here('preliminary_figures','settlement_counts_allregions.pdf'),width = 7,height = 6)
plot(sitesFreqAllRegions,type='n',xlab='BC/AD',ylab='site counts',xlim=c(-2100,250),axes=F)
polygon(x=c(BPtoBCAD(sitesFreqAllRegions$BP),BPtoBCAD(rev(sitesFreqAllRegions$BP))),y=c(sitesFreqAllRegions$Cnts,rep(0,nrow(sitesFreqAllRegions))),col='grey',border=NA)
axis(1,at=c(-2000,-1750,-1500,-1250,-1000,-750,-500,-250,1,250),labels=abs(c(-2000,-1750,-1500,-1250,-1000,-750,-500,-250,1,250)))
abline(v=c(-2000,-1750,-1500,-1250,-1000,-750,-500,-250,1,250),lty=3,col='white',lwd=2)
axis(2)
box()
dev.off()

## Per Region ----
regions=c('Hokkaido','Tohoku','Kanto','Chubu','Chugoku','Kansai','Shikoku','Kyushu')
cntByRegion = vector('list',length = length(regions))
names(cntByRegion) = regions
for (i in 1:length(regions))
{
  cntByRegion[[i]] = countSites(ageRange=timeRange,ageMatrix = ageMatrix,sites=subset(sites,Region==regions[i]),ageRef=ages)
}

pdf(file = here('preliminary_figures','settlement_counts_per_region.pdf'),width = 12,height = 6)
par(mfrow=c(2,4),mar=c(6,5,3,1))

for(i in 1:length(cntByRegion))
{
  tmp=cntByRegion[[i]]
  plot(sitesFreqAllRegions,type='n',xlab='BC/AD',ylab='site counts',xlim=c(-2100,250),ylim=range(tmp$Cnts),axes=F)
  polygon(x=c(BPtoBCAD(tmp$BP),BPtoBCAD(rev(tmp$BP))),y=c(tmp$Cnts,rep(0,nrow(tmp))),col='grey',border=NA)
  axis(1,at=c(-2000,-1750,-1500,-1250,-1000,-750,-500,-250,1,250),labels=abs(c(-2000,-1750,-1500,-1250,-1000,-750,-500,-250,1,250)))
  axis(2)
  #abline(v=c(-2000,-1750,-1500,-1250,-1000,-750,-500,-250,1,250),lty=3,col='white',lwd=2)
  box()
  legend('topleft',legend=c(regions[i]),bty='n',cex=1.4)
}
dev.off()

# Permutation Test ----
sites=subset(sites,!is.na(Region))
permRes=countSitePermTest(ageRange=timeRange,ageMatrix = ageMatrix,sites=sites,ageRef=ages,nsim=500,ncores=2,backsight=200,regions=regions)
save(permRes,file=here('R_images','permRes_settlement.RData'))

## Plot Results ---- 
## Needs conversion into BCAD
pdf(file = here('preliminary_figures','permTest_counts_perRegion.pdf'),width = 12,height = 6)
par(mfrow=c(2,4),mar=c(6,5,3,1))

for(i in 1:length(regions))
{
  plotPermResult(permRes,method=c('counts'),index=i,calendar='BCAD',main=regions[i])
}
dev.off()