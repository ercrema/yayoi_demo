# Load Data and Library ----
library(here)
library(coda)
load(here('results','icar_c14doubleRes500.RData'))
load(here('results','icar_c14doubleRes750.RData'))

# Convert growth rate into % growth
icar.c14double500[,grep('s1\\[',colnames(icar.c14double500))]  <- icar.c14double500[,grep('s1\\[',colnames(icar.c14double500))] * 100
icar.c14double500[,grep('s2\\[',colnames(icar.c14double500))]  <- icar.c14double500[,grep('s2\\[',colnames(icar.c14double500))] * 100
icar.c14double750[,grep('s1\\[',colnames(icar.c14double750))]  <- icar.c14double750[,grep('s1\\[',colnames(icar.c14double750))] * 100
icar.c14double750[,grep('s2\\[',colnames(icar.c14double750))]  <- icar.c14double750[,grep('s2\\[',colnames(icar.c14double750))] * 100

# Transform sigma into tau ----
icar.c14double500$tau1  <- 1/sqrt(icar.c14double500$sigma1)
icar.c14double500$tau2  <- 1/sqrt(icar.c14double500$sigma2)
icar.c14double750$tau1  <- 1/sqrt(icar.c14double750$sigma1)
icar.c14double750$tau2  <- 1/sqrt(icar.c14double750$sigma2)
icar.c14double500 <- icar.c14double500[,-c(17,18)]
icar.c14double750 <- icar.c14double750[,-c(17,18)]

# Parameter Names ----
params <- c(paste0('r_1_',as.roman(1:8)),paste0('r_2_',as.roman(1:8)),'tau_1','tau_2')
table.S1  <- table.S2  <-  data.frame(params)

# Mean Posteriors ----
table.S1$mean <- apply(icar.c14double500,2,mean)
table.S2$mean <- apply(icar.c14double500,2,mean)

# HPDI ----
table.S1$hpd <- apply(icar.c14double500,2,function(x){return(HPDinterval(as.mcmc(x),0.95))}) |> round(3) |> apply(2,paste,collapse=' ~ ')
table.S2$hpd <- apply(icar.c14double750,2,function(x){return(HPDinterval(as.mcmc(x),0.95))}) |> round(3) |> apply(2,paste,collapse=' ~ ')

# Rhat ----
table.S1$Rhat  <- rhats.c14double500[[1]][,1] |> round(4)
table.S2$Rhat  <- rhats.c14double750[[1]][,1] |> round(4)

# ESS ----
table.S1$ESS  <- ess.c14double500 |> round()
table.S2$ESS  <- ess.c14double750 |> round()

# Define Colnames and Export ----
colnames(table.S1)  <- colnames(table.S2)  <- c('Parameter','Mean Posterior','95% HPD Interval','Rhat','ESS')
write.csv(table.S1,here('tables','tableS1.csv'),row.names=FALSE)
write.csv(table.S2,here('tables','tableS2.csv'),row.names=FALSE)
