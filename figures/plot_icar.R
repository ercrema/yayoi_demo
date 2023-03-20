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
load(here('results','icar_finaljomon.RData'))
load(here('results','icar_yayoi.RData'))
load(here('results','icar_before.RData'))
load(here('results','icar_after.RData'))

# Extract summaries
finaljomon.m <- apply(icar.finaljomon,2,median) * 100
finaljomon.lo <- apply(icar.finaljomon,2,quantile,0.1) * 100
finaljomon.hi <- apply(icar.finaljomon,2,quantile,0.9) * 100

yayoi.m <- apply(icar.yayoi,2,median) * 100
yayoi.lo <- apply(icar.yayoi,2,quantile,0.1) * 100
yayoi.hi <- apply(icar.yayoi,2,quantile,0.9) * 100

before.m <- apply(icar.before,2,median) * 100
before.lo <- apply(icar.before,2,quantile,0.1) * 100
before.hi <- apply(icar.before,2,quantile,0.9) * 100

after.m <- apply(icar.after,2,median) * 100
after.lo <- apply(icar.after,2,quantile,0.1) * 100
after.hi <- apply(icar.after,2,quantile,0.9) * 100

limits  <- c(min(c(after.lo,finaljomon.lo)),max(c(after.hi,finaljomon.hi)))


# Load Base Map ----
win.sf <- ne_states(country = "japan",returnclass='sf') |> subset(!name_vi %in% c("Okinawa","Hokkaido"))

 

# Final Jomon Results ----
win.sf$pred_r_m = finaljomon.m[1:45]
win.sf$pred_r_lo = finaljomon.lo[1:45]
win.sf$pred_r_hi = finaljomon.hi[1:45]

g1<-ggplot(win.sf,aes(fill=pred_r_m))  + geom_sf() + scale_fill_gradient2(limits=limits)  + xlab("Longitude") + ylab("Latitude") + xlim(130,142) + ylim(29,42) + labs(fill='Percentage \n Growth Rate') + ggtitle('Median Posterior Growth Rate \n (2000-1000 BC)') + theme(legend.position = c(0.8, 0.2), plot.title=element_text(vjust=-15,hjust=0.1,size=20)) 
g2<-ggplot(win.sf,aes(fill=pred_r_lo))  + geom_sf() + scale_fill_gradient2(limits=limits) + xlab("Longitude") + ylab("Latitude") + xlim(130,142) + ylim(29,42) + ggtitle('10th Posterior \nPercentile') + theme(legend.position = "none",plot.title=element_text(vjust=-12,hjust=0.1))
g3<-ggplot(win.sf,aes(fill=pred_r_hi))  + geom_sf() + scale_fill_gradient2(limits=limits) + xlab("Longitude") + ylab("Latitude") + xlim(130,142) + ylim(29,42) + ggtitle('90th Posterior \nPercentile') + theme(legend.position = "none",plot.title=element_text(vjust=-12,hjust=0.1))

lay <- rbind(c(1,1,2),c(1,1,3))
pdf(here('figures','icarfinaljomon.pdf'),width=9,height=8)
grid.arrange(grobs=list(g1,g2,g3),layout_matrix=lay)
dev.off()


# Comparative Plot
icar.fj  <- icar.finaljomon[,1:45]
colnames(icar.fj)  <- unlist(lapply(strsplit(win.sf$name_en," "),function(x){x[1]}))
ii  <- order(win.sf$latitude)
icar.fj  <- icar.fj[,ii]
finaljomon  <- gather(icar.fj)
finaljomon$key  <- factor(finaljomon$key,levels=colnames(icar.fj))
finaljomon$value  <-  finaljomon$value * 100

finaljomon2  <- finaljomon
finaljomon2$value  <- finaljomon2$value + 0.01


g.finalj.bars <- ggplot(aes(y = key, x = value),data=finaljomon) + 
	     	 stat_interval(aes(y = key), .width = c(.5, .8, .95)) +
# 	         stat_interval(aes(y = key), .width = c(.5, .8, .95),data=finaljomon2, position=position_nudge(y = -0.5)) +
	         scale_color_brewer() +
	         ylab('')   +
	         xlab('Growth Rate (%)') +
#  	         ggtitle('Estimated Growth Rate (2000-1000 BC)') +
	         theme(legend.position=c(0.9,0.9)) +
	         geom_vline(xintercept=0, linetype='dashed')

g.finalj.map <-ggplot(win.sf,aes(fill=pred_r_m)) + 
	       geom_sf() + 
	       scale_fill_gradient2(limits=limits) +
	       xlab("Longitude") +
	       ylab("Latitude") +
	       xlim(130,142) +
	       ylim(29,42) +
	       labs(fill='Percentage \n Growth Rate') +
	       ggtitle('Median Posterior Growth Rate \n (2000-1000 BC)') +
	       theme(legend.position = c(0.8, 0.2), plot.title=element_text(size=20)) 
pdf(here('figures','icarfinaljomon2.pdf'),width=10,height=7)
grid.arrange(g.finalj.map,g.finalj.bars,ncol=2)
dev.off()

# Yayoi Results ----

win.sf$pred_r_m = yayoi.m[1:45]
win.sf$pred_r_lo = yayoi.lo[1:45]
win.sf$pred_r_hi = yayoi.hi[1:45]

g1<-ggplot(win.sf,aes(fill=pred_r_m))  + geom_sf() + scale_fill_gradient2(limits=limits)  + xlab("Longitude") + ylab("Latitude") + xlim(130,142) + ylim(29,42) + labs(fill='Percentage \n Growth Rate') + ggtitle('Median Posterior Growth Rate \n Yayoi period (1000BC - 250AD') + theme(legend.position = c(0.8, 0.2), plot.title=element_text(vjust=-15,hjust=0.1,size=20)) 
g2<-ggplot(win.sf,aes(fill=pred_r_lo))  + geom_sf() + scale_fill_gradient2(limits=limits) + xlab("Longitude") + ylab("Latitude") + xlim(130,142) + ylim(29,42) + ggtitle('10th Posterior \nPercentile') + theme(legend.position = "none",plot.title=element_text(vjust=-12,hjust=0.1))
g3<-ggplot(win.sf,aes(fill=pred_r_hi))  + geom_sf() + scale_fill_gradient2(limits=limits) + xlab("Longitude") + ylab("Latitude") + xlim(130,142) + ylim(29,42) + ggtitle('90th Posterior \nPercentile') + theme(legend.position = "none",plot.title=element_text(vjust=-12,hjust=0.1))

lay <- rbind(c(1,1,2),c(1,1,3))
pdf(here('figures','icaryayoi.pdf'),width=9,height=8)
grid.arrange(grobs=list(g1,g2,g3),layout_matrix=lay)
dev.off()


# Comparative Plot
icar.y  <- icar.yayoi[,1:45]
colnames(icar.y)  <- unlist(lapply(strsplit(win.sf$name_en," "),function(x){x[1]}))
ii  <- order(win.sf$latitude)
icar.y  <- icar.y[,ii]
yayoi  <- gather(icar.y)
yayoi$key  <- factor(yayoi$key,levels=colnames(icar.y))
yayoi$value  <-  yayoi$value * 100

g.yayoi.bars <- ggplot(aes(y = key, x = value),data=yayoi) + 
	     	 stat_interval(aes(y = key), .width = c(.5, .8, .95)) +
# 	         stat_interval(aes(y = key), .width = c(.5, .8, .95),data=yayoi2, position=position_nudge(y = -0.5)) +
	         scale_color_brewer() +
	         ylab('')   +
	         xlab('Growth Rate (%)') +
#  	         ggtitle('Estimated Growth Rate (2000-1000 BC)') +
	         theme(legend.position=c(0.9,0.5)) +
	         geom_vline(xintercept=0, linetype='dashed')

g.yayoi.map <-ggplot(win.sf,aes(fill=pred_r_m)) + 
	       geom_sf() + 
	       scale_fill_gradient2(limits=limits) +
	       xlab("Longitude") +
	       ylab("Latitude") +
	       xlim(130,142) +
	       ylim(29,42) +
	       labs(fill='Percentage \n Growth Rate') +
	       ggtitle('Median Posterior Growth Rate \n Yayoi Period (1000BC - 250AD)') +
	       theme(legend.position = c(0.8, 0.2), plot.title=element_text(size=20)) 
pdf(here('figures','icaryayoi2.pdf'),width=10,height=7)
grid.arrange(g.yayoi.map,g.yayoi.bars,ncol=2)
dev.off()


# Before Results ----

win.sf$pred_r_m = before.m[1:45]
win.sf$pred_r_lo = before.lo[1:45]
win.sf$pred_r_hi = before.hi[1:45]

g1<-ggplot(win.sf,aes(fill=pred_r_m))  + geom_sf() + scale_fill_gradient2(limits=limits)  + xlab("Longitude") + ylab("Latitude") + xlim(130,142) + ylim(29,42) + labs(fill='Percentage \n Growth Rate') + ggtitle('Median Posterior Growth Rate \n 1000yrs before rice arrival') + theme(legend.position = c(0.8, 0.2), plot.title=element_text(vjust=-15,hjust=0.1,size=20)) 
g2<-ggplot(win.sf,aes(fill=pred_r_lo))  + geom_sf() + scale_fill_gradient2(limits=limits) + xlab("Longitude") + ylab("Latitude") + xlim(130,142) + ylim(29,42) + ggtitle('10th Posterior \nPercentile') + theme(legend.position = "none",plot.title=element_text(vjust=-12,hjust=0.1))
g3<-ggplot(win.sf,aes(fill=pred_r_hi))  + geom_sf() + scale_fill_gradient2(limits=limits) + xlab("Longitude") + ylab("Latitude") + xlim(130,142) + ylim(29,42) + ggtitle('90th Posterior \nPercentile') + theme(legend.position = "none",plot.title=element_text(vjust=-12,hjust=0.1))

lay <- rbind(c(1,1,2),c(1,1,3))
pdf(here('figures','icarbefore.pdf'),width=9,height=8)
grid.arrange(grobs=list(g1,g2,g3),layout_matrix=lay)
dev.off()


# Comparative Plot
icar.bf  <- icar.before[,1:45]
colnames(icar.bf)  <- unlist(lapply(strsplit(win.sf$name_en," "),function(x){x[1]}))
ii  <- order(win.sf$latitude)
icar.bf  <- icar.bf[,ii]
before  <- gather(icar.bf)
before$key  <- factor(before$key,levels=colnames(icar.bf))
before$value  <-  before$value * 100

g.before.bars <- ggplot(aes(y = key, x = value),data=before) + 
	     	 stat_interval(aes(y = key), .width = c(.5, .8, .95)) +
# 	         stat_interval(aes(y = key), .width = c(.5, .8, .95),data=before2, position=position_nudge(y = -0.5)) +
	         scale_color_brewer() +
	         ylab('')   +
	         xlab('Growth Rate (%)') +
#  	         ggtitle('Estimated Growth Rate (2000-1000 BC)') +
	         theme(legend.position=c(0.9,0.5)) +
	         geom_vline(xintercept=0, linetype='dashed')

g.before.map <-ggplot(win.sf,aes(fill=pred_r_m)) + 
	       geom_sf() + 
	       scale_fill_gradient2(limits=limits) +
	       xlab("Longitude") +
	       ylab("Latitude") +
	       xlim(130,142) +
	       ylim(29,42) +
	       labs(fill='Percentage \n Growth Rate') +
	       ggtitle('Median Posterior Growth Rate \n 1000 years before rice arrival') +
	       theme(legend.position = c(0.8, 0.2), plot.title=element_text(size=20)) 
pdf(here('figures','icarbefore2.pdf'),width=10,height=7)
grid.arrange(g.before.map,g.before.bars,ncol=2)
dev.off()

# After Results ----
win.sf$pred_r_m = after.m[1:45]
win.sf$pred_r_lo = after.lo[1:45]
win.sf$pred_r_hi = after.hi[1:45]

g1<-ggplot(win.sf,aes(fill=pred_r_m))  + geom_sf() + scale_fill_gradient2(limits=limits)  + xlab("Longitude") + ylab("Latitude") + xlim(130,142) + ylim(29,42) + labs(fill='Percentage \n Growth Rate') + ggtitle('Median Posterior Growth Rate \n 1000yrs after rice arrival') + theme(legend.position = c(0.8, 0.2), plot.title=element_text(vjust=-15,hjust=0.1,size=20)) 
g2<-ggplot(win.sf,aes(fill=pred_r_lo))  + geom_sf() + scale_fill_gradient2(limits=limits) + xlab("Longitude") + ylab("Latitude") + xlim(130,142) + ylim(29,42) + ggtitle('10th Posterior \nPercentile') + theme(legend.position = "none",plot.title=element_text(vjust=-12,hjust=0.1))
g3<-ggplot(win.sf,aes(fill=pred_r_hi))  + geom_sf() + scale_fill_gradient2(limits=limits) + xlab("Longitude") + ylab("Latitude") + xlim(130,142) + ylim(29,42) + ggtitle('90th Posterior \nPercentile') + theme(legend.position = "none",plot.title=element_text(vjust=-12,hjust=0.1))

lay <- rbind(c(1,1,2),c(1,1,3))
pdf(here('figures','icarafter.pdf'),width=9,height=8)
grid.arrange(grobs=list(g1,g2,g3),layout_matrix=lay)
dev.off()


# Comparative Plot
icar.af  <- icar.after[,1:45]
colnames(icar.af)  <- unlist(lapply(strsplit(win.sf$name_en," "),function(x){x[1]}))
ii  <- order(win.sf$latitude)
icar.af  <- icar.af[,ii]
after  <- gather(icar.af)
after$key  <- factor(after$key,levels=colnames(icar.af))
after$value  <-  after$value * 100

g.after.bars <- ggplot(aes(y = key, x = value),data=after) + 
	     	 stat_interval(aes(y = key), .width = c(.5, .8, .95)) +
# 	         stat_interval(aes(y = key), .width = c(.5, .8, .95),data=after2, position=position_nudge(y = -0.5)) +
	         scale_color_brewer() +
	         ylab('')   +
	         xlab('Growth Rate (%)') +
#  	         ggtitle('Estimated Growth Rate (2000-1000 BC)') +
	         theme(legend.position=c(0.9,0.5)) +
	         geom_vline(xintercept=0, linetype='dashed')

g.after.map <-ggplot(win.sf,aes(fill=pred_r_m)) + 
	       geom_sf() + 
	       scale_fill_gradient2(limits=limits) +
	       xlab("Longitude") +
	       ylab("Latitude") +
	       xlim(130,142) +
	       ylim(29,42) +
	       labs(fill='Percentage \n Growth Rate') +
	       ggtitle('Median Posterior Growth Rate \n 1000 years from arrival Farming') +
	       theme(legend.position = c(0.8, 0.2), plot.title=element_text(size=20)) 
pdf(here('figures','icarafter2.pdf'),width=10,height=7)
grid.arrange(g.after.map,g.after.bars,ncol=2)
dev.off()

# After Results ----
win.sf$pred_r_m = after.m[1:45]
win.sf$pred_r_lo = after.lo[1:45]
win.sf$pred_r_hi = after.hi[1:45]

g1<-ggplot(win.sf,aes(fill=pred_r_m))  + geom_sf() + scale_fill_gradient2(limits=limits)  + xlab("Longitude") + ylab("Latitude") + xlim(130,142) + ylim(29,42) + labs(fill='Percentage \n Growth Rate') + ggtitle('Median Posterior Growth Rate \n 1000yrs after rice arrival') + theme(legend.position = c(0.8, 0.2), plot.title=element_text(vjust=-15,hjust=0.1,size=20)) 
g2<-ggplot(win.sf,aes(fill=pred_r_lo))  + geom_sf() + scale_fill_gradient2(limits=limits) + xlab("Longitude") + ylab("Latitude") + xlim(130,142) + ylim(29,42) + ggtitle('10th Posterior \nPercentile') + theme(legend.position = "none",plot.title=element_text(vjust=-12,hjust=0.1))
g3<-ggplot(win.sf,aes(fill=pred_r_hi))  + geom_sf() + scale_fill_gradient2(limits=limits) + xlab("Longitude") + ylab("Latitude") + xlim(130,142) + ylim(29,42) + ggtitle('90th Posterior \nPercentile') + theme(legend.position = "none",plot.title=element_text(vjust=-12,hjust=0.1))

lay <- rbind(c(1,1,2),c(1,1,3))
pdf(here('figures','icarafter.pdf'),width=9,height=8)
grid.arrange(grobs=list(g1,g2,g3),layout_matrix=lay)
dev.off()


# Comparative Plot
icar.af  <- icar.after[,1:45]
colnames(icar.af)  <- unlist(lapply(strsplit(win.sf$name_en," "),function(x){x[1]}))
ii  <- order(win.sf$latitude)
icar.af  <- icar.af[,ii]
after  <- gather(icar.af)
after$key  <- factor(after$key,levels=colnames(icar.af))
after$value  <-  after$value * 100

after2  <- after
after2$value  <- after2$value + 0.01


g.after.bars <- ggplot(aes(y = key, x = value),data=after) + 
	     	 stat_interval(aes(y = key), .width = c(.5, .8, .95)) +
# 	         stat_interval(aes(y = key), .width = c(.5, .8, .95),data=after2, position=position_nudge(y = -0.5)) +
	         scale_color_brewer() +
	         ylab('')   +
	         xlab('Growth Rate (%)') +
#  	         ggtitle('Estimated Growth Rate (2000-1000 BC)') +
	         theme(legend.position=c(0.9,0.5)) +
	         geom_vline(xintercept=0, linetype='dashed')

g.after.map <-ggplot(win.sf,aes(fill=pred_r_m)) + 
	       geom_sf() + 
	       scale_fill_gradient2(limits=limits) +
	       xlab("Longitude") +
	       ylab("Latitude") +
	       xlim(130,142) +
	       ylim(29,42) +
	       labs(fill='Percentage \n Growth Rate') +
	       ggtitle('Median Posterior Growth Rate \n 1000 years from arrival Farming') +
	       theme(legend.position = c(0.8, 0.2), plot.title=element_text(size=20)) 
pdf(here('figures','icarafter2.pdf'),width=10,height=7)
grid.arrange(g.after.map,g.after.bars,ncol=2)
dev.off()

# After Results ----
win.sf$pred_r_m = after.m[1:45]
win.sf$pred_r_lo = after.lo[1:45]
win.sf$pred_r_hi = after.hi[1:45]

g1<-ggplot(win.sf,aes(fill=pred_r_m))  + geom_sf() + scale_fill_gradient2(limits=limits)  + xlab("Longitude") + ylab("Latitude") + xlim(130,142) + ylim(29,42) + labs(fill='Percentage \n Growth Rate') + ggtitle('Median Posterior Growth Rate \n 1000yrs after rice arrival') + theme(legend.position = c(0.8, 0.2), plot.title=element_text(vjust=-15,hjust=0.1,size=20)) 
g2<-ggplot(win.sf,aes(fill=pred_r_lo))  + geom_sf() + scale_fill_gradient2(limits=limits) + xlab("Longitude") + ylab("Latitude") + xlim(130,142) + ylim(29,42) + ggtitle('10th Posterior \nPercentile') + theme(legend.position = "none",plot.title=element_text(vjust=-12,hjust=0.1))
g3<-ggplot(win.sf,aes(fill=pred_r_hi))  + geom_sf() + scale_fill_gradient2(limits=limits) + xlab("Longitude") + ylab("Latitude") + xlim(130,142) + ylim(29,42) + ggtitle('90th Posterior \nPercentile') + theme(legend.position = "none",plot.title=element_text(vjust=-12,hjust=0.1))

lay <- rbind(c(1,1,2),c(1,1,3))
pdf(here('figures','icarafter.pdf'),width=9,height=8)
grid.arrange(grobs=list(g1,g2,g3),layout_matrix=lay)
dev.off()


# Comparative Plot
icar.af  <- icar.after[,1:45]
colnames(icar.af)  <- unlist(lapply(strsplit(win.sf$name_en," "),function(x){x[1]}))
ii  <- order(win.sf$latitude)
icar.af  <- icar.af[,ii]
after  <- gather(icar.af)
after$key  <- factor(after$key,levels=colnames(icar.af))
after$value  <-  after$value * 100

after2  <- after
after2$value  <- after2$value + 0.01


g.after.bars <- ggplot(aes(y = key, x = value),data=after) + 
	     	 stat_interval(aes(y = key), .width = c(.5, .8, .95)) +
# 	         stat_interval(aes(y = key), .width = c(.5, .8, .95),data=after2, position=position_nudge(y = -0.5)) +
	         scale_color_brewer() +
	         ylab('')   +
	         xlab('Growth Rate (%)') +
#  	         ggtitle('Estimated Growth Rate (2000-1000 BC)') +
	         theme(legend.position=c(0.9,0.5)) +
	         geom_vline(xintercept=0, linetype='dashed')

g.after.map <-ggplot(win.sf,aes(fill=pred_r_m)) + 
	       geom_sf() + 
	       scale_fill_gradient2(limits=limits) +
	       xlab("Longitude") +
	       ylab("Latitude") +
	       xlim(130,142) +
	       ylim(29,42) +
	       labs(fill='Percentage \n Growth Rate') +
	       ggtitle('Median Posterior Growth Rate \n 1000 years from arrival Farming') +
	       theme(legend.position = c(0.8, 0.2), plot.title=element_text(size=20)) 
pdf(here('figures','icarafter2.pdf'),width=10,height=7)
grid.arrange(g.after.map,g.after.bars,ncol=2)
dev.off()

