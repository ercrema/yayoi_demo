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
regionPrefMatch$Reg[which(regionPrefMatch$Reg %in% c('Chugoku','Shikoku'))] = 'Chugoku + Shikoku'
regionPrefMatch$Reg[which(regionPrefMatch$Pref %in% c('Niigata','Fukushima','Miyagi','Yamagata'))] = 'S. Tohoku + Niigata'
regionPrefMatch$Reg[which(regionPrefMatch$Reg == 'Tohoku')] = 'N. Tohoku'
regionPrefMatch$ID <- NA
regionPrefMatch$ID[which(regionPrefMatch$Reg == 'Kyushu')] = 1
regionPrefMatch$ID[which(regionPrefMatch$Reg == 'Chugoku + Shikoku')] = 2
regionPrefMatch$ID[which(regionPrefMatch$Reg == 'Kansai')] = 3
regionPrefMatch$ID[which(regionPrefMatch$Reg == 'Chubu')] = 4
regionPrefMatch$ID[which(regionPrefMatch$Reg == 'Kanto')] = 5
regionPrefMatch$ID[which(regionPrefMatch$Reg == 'S. Tohoku + Niigata')] = 6
regionPrefMatch$ID[which(regionPrefMatch$Reg == 'N. Tohoku')] = 7
regionPrefMatch$ID[which(regionPrefMatch$Reg == 'Hokkaido')] = 8

win@data <- left_join(win@data,regionPrefMatch,by=c('name'='Pref'))
win2 <- gUnaryUnion(win,id=win@data$ID)
win.sf <- as(win,'sf')
win2.sf <- as(win2,'sf')

# Read ICAR Results
load(here('R_images','icarres.RData'))
samples.m <- apply(icar.samples[[1]],2,median)
samples.lo <- apply(icar.samples[[1]],2,quantile,0.025)
samples.hi <- apply(icar.samples[[1]],2,quantile,0.975)

win.sf <- as(win,'sf')
win.sf$pred_r_m = samples.m[1:45] *100
win.sf$pred_r_lo = samples.lo[1:45] *100
win.sf$pred_r_hi = samples.hi[1:45] *100
library(gridExtra)
pdf(file=here('figures','icar_demo_res.pdf'),width=10,height=8)
g1<-ggplot()  + geom_sf(win.sf, mapping=aes(fill=pred_r_m),size=.05) + scale_fill_gradient2(limits=c(-0.2,0.41)) + ggtitle('Posterior Median') + xlab("Longitude") + ylab("Latitude") + xlim(130,142) + ylim(29,42)  + geom_sf(data=win2.sf,size=0.8,fill=NA) + labs(fill='Annual Growth Rate (%)') + theme(legend.position=c(0.2,0.8))
g2<-ggplot()  + geom_sf(win.sf, mapping=aes(fill=pred_r_lo),size=.05) + scale_fill_gradient2(limits=c(-0.2,0.41)) + ggtitle('Posterior 2.5th Percentile') + xlab("Longitude") + ylab("Latitude") + xlim(130,142) + ylim(29,42)  + geom_sf(data=win2.sf,size=0.8,fill=NA) + labs(fill='') + theme(axis.text=element_blank(),legend.position='none')
g3<-ggplot()  + geom_sf(win.sf, mapping=aes(fill=pred_r_hi),size=.05) + scale_fill_gradient2(limits=c(-0.2,0.41)) + ggtitle('Predicted 97.5th Percentile') + xlab("Longitude") + ylab("Latitude") + xlim(130,142) + ylim(29,42)  + geom_sf(data=win2.sf,size=0.8,fill=NA) + labs(fill='') + theme(axis.text=element_blank(),legend.position='none')
lay <- rbind(c(1,1,2),c(1,1,3))
grid.arrange(g1,g2,g3,layout_matrix=lay)
dev.off()
