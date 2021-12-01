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

# Time window of analyses ----
a <- BCADtoBP(-150)
b <- BCADtoBP(250)


# Load and prepare 14C and Spatial Data ----
# Load Data
load(here('data','c14data.RData'))
c14db = subset(c14db,C14Age<3500 & C14Age> 1500 &  Material == 'Terrestrial' & !PrefectureNameEn %in% c('Hokkaido','Okinawa')) |> select(C14Age,C14Error,SiteID,Prefecture=PrefectureNameEn,Longitude,Latitude)

# Consider dates with probability mass above 0.5 within window of analyses
calibrated.dates <- calibrate(c14db$C14Age,c14db$C14Error)
i <- which.CalDates(calibrated.dates,BP<a&BP>b,p=0.5)
c14db <- c14db[i,]
calibrated.dates <- calibrated.dates[i]

# Site Level Thinning
bin50  <- binPrep(sites = c14db$SiteID,ages=calibrated.dates,h=50)
i  <- thinDates(ages=c14db$C14Age,  errors=c14db$C14Error, bins=bin50, size=1, thresh=1,seed=123,method='splitsample')
c14db  <-  c14db[i,]
calibrated.dates <- calibrated.dates[i]

# Define initialisation values for theta
theta.init <- medCal(calibrated.dates)
theta.init[which(theta.init>=a)] = a - 1
theta.init[which(theta.init<=b)] = b + 1

# Define total number of dates
N <- nrow(c14db)

# Define Spatial Structure
win <- ne_states(country = "japan") |> subset(!name_vi %in% c("Okinawa","Hokkaido"))
# sites <- select(c14db,SiteID,Longitude,Latitude)
# coordinates(sites) <- c('Longitude','Latitude')
# proj4string(sites) <- proj4string(win)
# sites.sf <- as(sites,'sf')

# Create Prefecture ID and add to c14db
Npref <- length(win)
win@data$ID = 1:Npref
c14db$PrefID = match(c14db$Prefecture,win@data$name)
sort(table(c14db$PrefID))

# Test Plot
# win.sf = as(win,'sf')
# ggplot(data=win.sf) +  geom_sf() + xlab("Longitude") + ylab("Latitude") + xlim(130,142) + ylim(29,42) + geom_sf(data = sites.sf,shape=20,size=0.5)

# Define Neighbourhood Structure
W_nb <- poly2nb(win, row.names =  win@data$ID)
# Determine neighborhood/adjacency information needed for neighborhood-based CAR model
nbInfo <- nb2WB(W_nb)

# Define Data
d <- list(cra=c14db$C14Age,cra_error=c14db$C14Error)


# Define Constants
data(intcal20)
constants  <-  list(Npref=Npref,N=N,adj=nbInfo$adj,weights=nbInfo$weights,num=nbInfo$num,L=length(nbInfo$adj),calBP=intcal20$CalBP,C14BP=intcal20$C14Age,C14err=intcal20$C14Age.sigma)
constants$a <- a
constants$b <- b
constants$id.pref <-c14db$PrefID



# Core Analysis ----

runFun <- function(seed, d, constants, theta.init, nburnin, niter, thin)
{
	library(nimbleCarbon)
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

	## Setup Init
	set.seed(seed)
	inits <- list(sigma = runif(1,0,100), s = rnorm(constants$Npref,sd=0.001),theta=theta.init)

	#MCMC
	model <- nimbleModel(icarmodel, constants = constants, data = d, inits = inits)
	cModel <- compileNimble(model)
	conf <- configureMCMC(model, monitors = c('s','sigma'))
	MCMC <- buildMCMC(conf)
	cMCMC <- compileNimble(MCMC, project = cModel)
	samples <- runMCMC(cMCMC, niter = niter, thin=thin,nburnin = nburnin,samplesAsCodaMCMC = T, setSeed = seed)
	return(samples)
}

# Run MCMC ----
ncores  <-  3
cl <- makeCluster(ncores)
seeds  <-  c(123,456,789)
niter  <- 250000
nburnin  <- 125000
thin  <- 5

chain_output  <- parLapply(cl = cl, X = seeds, fun = runFun, d = d,constants = constants, theta = theta.init, niter = niter, nburnin = nburnin,thin = thin)
stopCluster(cl)
icar.samples=coda::mcmc.list(chain_output)
coda::gelman.diag(icar.samples)
save(icar.samples,file=here('R_images','icarres_150to250.RData'))

# samples.m <- apply(icar.samples[[1]],2,median)
# samples.lo <- apply(icar.samples[[1]],2,quantile,0.025)
# samples.hi <- apply(icar.samples[[1]],2,quantile,0.975)
# 
# win.sf <- as(win,'sf')
# win.sf$pred_r_m = samples.m[1:45]
# win.sf$pred_r_lo = samples.lo[1:45]
# win.sf$pred_r_hi = samples.hi[1:45]
# library(gridExtra)
# g1<-ggplot(win.sf,aes(fill=pred_r_m))  + geom_sf() + scale_fill_gradient2() + ggtitle('Growth Rates') + xlab("Longitude") + ylab("Latitude") + xlim(130,142) + ylim(29,42)
# g2<-ggplot(win.sf,aes(fill=pred_r_lo))  + geom_sf() + scale_fill_gradient2() + ggtitle('Growth Rates') + xlab("Longitude") + ylab("Latitude") + xlim(130,142) + ylim(29,42)
# g3<-ggplot(win.sf,aes(fill=pred_r_hi))  + geom_sf() + scale_fill_gradient2() + ggtitle('Growth Rates') + xlab("Longitude") + ylab("Latitude") + xlim(130,142) + ylim(29,42)
# 
# lay <- rbind(c(1,1,2),c(1,1,3))
# grid.arrange(g1,g2,g3,layout_matrix=lay)

