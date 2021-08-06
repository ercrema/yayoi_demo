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

# General Parameters ----
a <- 3000
b <- 1700


# Load and prepare 14C and Spatial Data ----
# C14 Data between 3000 and 1700 calBP
load(here('data','c14data.RData'))
c14db <- subset(c14db,retain == TRUE & anthropic==TRUE & Material=='Terrestrial' & !PrefectureNameEn %in% c('Hokkaido','Okinawa'))
calibrated <- calibrate(c14db$CRA,c14db$CRAError)
ii <- which.CalDates(calibrated,BP<a&BP>b,p=0.5)
# possibly random thinning too? #reduction to anthrophic drastically reduces sample sizes...
c14db <- c14db[ii,]
calibrated <- calibrated[ii]
theta.init <- medCal(calibrated)
theta.init[which(theta.init>=a)] = a - 1
theta.init[which(theta.init<=b)] = b + 1
# table(c14db$PrefectureNameEn)
N <- nrow(c14db)





# Spatial Data
win <- ne_states(country = "japan") |> subset(!name_vi %in% c("Okinawa","Hokkaido"))
sites <- select(c14db,SiteID,Longitude,Latitude)
coordinates(sites) <- c('Longitude','Latitude')
proj4string(sites) <- proj4string(win)
sites.sf <- as(sites,'sf')

# Create Prefecture ID and add to c14db
Npref <- length(win)
win@data$ID = 1:Npref
c14db$PrefID = match(c14db$PrefectureNameEn,win@data$name)

# Test Plot
# win.sf = as(win,'sf')
# ggplot(data=win.sf) +  geom_sf() + xlab("Longitude") + ylab("Latitude") + xlim(130,142) + ylim(29,42) + geom_sf(data = sites.sf,shape=20,size=0.5)

# Define Neighbourhood Structure
W_nb <- poly2nb(win, row.names =  win@data$ID)
# Determine neighborhood/adjacency information needed for neighborhood-based CAR model
nbInfo <- nb2WB(W_nb)


# Core Analysis ----
## ICAR  Model
icarmodel <- nimbleCode({
	for (i in 1:N)
	{
		theta[i] ~ dExponentialGrowth(a=a,b=b,r=s[id.pref[i]])
		c14age[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
		sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
		sigmaDate[i] <- (cra_error[i]^2+sigmaCurve[i]^2)^(1/2);
		cra[i] ~ dnorm(mean=c14age[i],sd=sigmaDate[i]);
	}

	# ICAR Model Prior
	s[1:Npref] ~ dcar_normal(adj[1:L], weights[1:L], num[1:Npref], tau, zero_mean =0)
	tau <- 1/sigma^2
	sigma ~ dunif(0,100)
})

## Setup Constants
data(intcal20)
constants  <-  list(Npref=Npref,N=N,adj=nbInfo$adj,weights=nbInfo$weights,num=nbInfo$num,L=length(nbInfo$adj),calBP=intcal20$CalBP,C14BP=intcal20$C14Age,C14err=intcal20$C14Age.sigma)
constants$a <- a
constants$b <- b
constants$id.pref <-c14db$PrefID


## Setup Data
data <- list(cra=c14db$CRA,cra_error=c14db$CRAError) 

## Setup Init
inits <- list(sigma = 0.001, s = rnorm(Npref,sd=0.001),theta=theta.init)

## Run
model <- nimbleModel(icarmodel, constants = constants, data = data, inits = inits)
cModel <- compileNimble(model)
conf <- configureMCMC(model, monitors = c('s','sigma'))
conf$printSamplers()
MCMC <- buildMCMC(conf)
cMCMC <- compileNimble(MCMC, project = cModel)
samples <- runMCMC(cMCMC, niter = 10000, nburnin = 5000, thin = 1)
samples.mu <- apply(samples,2,mean)
win.sf$pred_r = samples.mu[1:45]
ggplot(win.sf,aes(fill=pred_r))  + geom_sf() + scale_fill_gradient2() + ggtitle('Growth Rates') + xlab("Longitude") + ylab("Latitude") + xlim(130,142) + ylim(29,42)



