library(ggplot2)
library(dplyr)
library(gridExtra)
library(nimbleCarbon)
library(rcarbon)
library(rnaturalearth)
library(sf)
library(maptools)
library(rgeos)
library(spdep)
library(here)
library(parallel)

# Read Spatial Data ----
win <- ne_states(country = "japan") |> subset(!name_vi %in% c("Okinawa","Hokkaido"))

load(here('data','japanc14db_v04.1(210715).RData'))
regionPrefMatch  <- select(c14db,Pref=PrefectureNameEn,Reg=Region) |> unique()
regionPrefMatch <- rbind.data.frame(regionPrefMatch,data.frame(Pref='Yamaguchi',Reg='Chugoku'))

win@data <- left_join(win@data,regionPrefMatch,by=c('name'='Pref'))
win2 <- gUnaryUnion(win,id=win@data$ID)
win.sf <- as(win,'sf')
win2.sf <- as(win2,'sf')

# Read ICAR Results
win.sf <- as(win,'sf')
load(here('R_images','icarres_950to550.RData'))
coda::gelman.diag(icar.samples)
samples1.m <- apply(icar.samples[[1]],2,mean,na.rm=T)
samples1.lo <- apply(icar.samples[[1]],2,quantile,0.025,na.rm=T)
samples1.hi <- apply(icar.samples[[1]],2,quantile,0.975,na.rm=T)
win.sf$pred_r_m1 = samples1.m[1:45] *100
win.sf$pred_r_lo1 = samples1.lo[1:45] *100
win.sf$pred_r_hi1 = samples1.hi[1:45] *100

load(here('R_images','icarres_550to150.RData'))
coda::gelman.diag(icar.samples)
samples2.m <- apply(icar.samples[[1]],2,mean,na.rm=T)
samples2.lo <- apply(icar.samples[[1]],2,quantile,0.025,na.rm=T)
samples2.hi <- apply(icar.samples[[1]],2,quantile,0.975,na.rm=T)
win.sf$pred_r_m2 = samples2.m[1:45] *100
win.sf$pred_r_lo2 = samples2.lo[1:45] *100
win.sf$pred_r_hi2 = samples2.hi[1:45] *100

load(here('R_images','icarres_150to250.RData'))
coda::gelman.diag(icar.samples)
samples3.m <- apply(icar.samples[[1]],2,mean,na.rm=T)
samples3.lo <- apply(icar.samples[[1]],2,quantile,0.025,na.rm=T)
samples3.hi <- apply(icar.samples[[1]],2,quantile,0.975,na.rm=T)
win.sf$pred_r_m3 = samples3.m[1:45] *100
win.sf$pred_r_lo3 = samples3.lo[1:45] *100
win.sf$pred_r_hi3 = samples3.hi[1:45] *100


library(gridExtra)
pdf(file=here('figures','icar_demo_res.pdf'),width=10,height=8)

g1<-ggplot()  + geom_sf(win.sf, mapping=aes(fill=pred_r_m1),size=.05) + scale_fill_gradient2(limits=c(-0.2,0.6)) + ggtitle('Posterior Median') + xlab("Longitude") + ylab("Latitude") + xlim(130,142) + ylim(29,42)  +  labs(fill='Annual Growth Rate (%)') + theme(legend.position=c(0.2,0.8))
g2<-ggplot()  + geom_sf(win.sf, mapping=aes(fill=pred_r_m2),size=.05) + scale_fill_gradient2(limits=c(-0.2,0.6)) + ggtitle('Posterior Median') + xlab("Longitude") + ylab("Latitude") + xlim(130,142) + ylim(29,42)  +  labs(fill='Annual Growth Rate (%)') + theme(legend.position=c(0.2,0.8))
g3<-ggplot()  + geom_sf(win.sf, mapping=aes(fill=pred_r_m3),size=.05) + scale_fill_gradient2(limits=c(-0.2,0.6)) + ggtitle('Posterior Median') + xlab("Longitude") + ylab("Latitude") + xlim(130,142) + ylim(29,42)  +  labs(fill='Annual Growth Rate (%)') + theme(legend.position=c(0.2,0.8))
