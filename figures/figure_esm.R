library(here)
library(RColorBrewer)


# Figure S1 #Scatter plot Site vs Pop Estimate Haji ----
comp.est  <- read.csv(here('results','pop_estimate_compare.csv'))
pdf(width=8,height=6,file=here('figures','figure_S1.pdf'))
plot(comp.est$Haji.pop,comp.est$Haji.sites,pch=20,xlab='Population Estimate',ylab='Number of Sites',main='Comparison Site Number of Historical Population Estimate for Haji period',ylim=c(0,40500),xlim=c(270000,1230000))
text(comp.est$Haji.pop,comp.est$Haji.site+1000,labels=comp.est$Region)
dev.off()

# Figure S2 #Comparison estimates ----
comp.est  <- comp.est[match(c('Kyushu','Shikoku','Chugoku','Kinki','Tokai','Chubu','Hokuriku','Kanto','Tohoku'),comp.est$Region),]
cols  <- brewer.pal(4,'Set2')
pdf(width=8,height=6,file=here('figures','figure_S2.pdf'))
barplot(t(as.matrix(comp.est[,c('Yayoi.pop84','Yayoi.pop.est1','Yayoi.pop.est2','Yayoi.pop.est3')])),beside=TRUE,names.arg=comp.est$Region,col=cols,ylab='Population Size',cex.names=0.8)
legend('topright',legend=c('Koyama and Sugito 1984','Estimate 1','Estimate 2','Estimate 3'),fill=cols,bty='n')
dev.off()



