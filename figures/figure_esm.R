library(here)
library(nimbleCarbon)
library(RColorBrewer)
library(latex2exp)

# Figure S1 #Double Exponential Model ----
r1=c(0.01,0.001,-0.005)
r2=c(0.01,0.001,-0.005)
params = expand.grid(r1=r1,r2=r2)

pdf(width=8,height=8,file=here('figures','figure_S1.pdf'))
par(mfrow=c(3,3))
for (i in 1:nrow(params))
{
	modelPlot(dDoubleExponentialGrowth,a=3000,b=2000,params=c(mu=2500,r1=params$r1[i],r2=params$r2[i]),col=1,alpha=1)
	abline(v=2500,lty=2)
	title(paste('r_1=',params$r1[i],';r_2=',params$r2[i],sep=''))
}
dev.off()

# Figure S2 #Scatter plot Site vs Pop Estimate Haji ----
comp.est  <- read.csv(here('results','pop_estimate_compare.csv'))
pdf(width=8,height=6,file=here('figures','figure_S2.pdf'))
plot(comp.est$Haji.pop,comp.est$Haji.sites,pch=20,xlab='Population Estimate',ylab='Number of Sites',main='Comparison Site Number of Historical Population Estimate for Haji period',ylim=c(0,40500),xlim=c(270000,1230000))
text(comp.est$Haji.pop,comp.est$Haji.site+1000,labels=comp.est$Region)
dev.off()

# Figure S3 #Comparison estimates ----
comp.est  <- comp.est[match(c('Kyushu','Shikoku','Chugoku','Kinki','Tokai','Chubu','Hokuriku','Kanto','Tohoku'),comp.est$Region),]
cols  <- brewer.pal(4,'Set2')
pdf(width=8,height=6,file=here('figures','figure_S3.pdf'))
barplot(t(as.matrix(comp.est[,c('Yayoi.pop84','Yayoi.pop.est1','Yayoi.pop.est2','Yayoi.pop.est3')])),beside=TRUE,names.arg=comp.est$Region,col=cols,ylab='Population Size',cex.names=0.8)
legend('topright',legend=c('Koyama and Sugito 1984','Estimate 1','Estimate 2','Estimate 3'),fill=cols,bty='n')
dev.off()



