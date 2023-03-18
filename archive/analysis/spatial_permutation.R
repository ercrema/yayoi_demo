# Load Library and Data ----
library(rcarbon)
library(here)
load(here('data','c14data.RData'))

# Preparation ----
BPinterval  <-  BCADtoBP(c(-1500,750))
c14db  <- subset(c14db,CRA<=BPinterval[1] & CRA >=BPinterval[2] & MaterialCode1 == 'T')

# Define Sites
sites  <- unique(data.frame(id=c14db$SiteID,lat=c14db$Latitude,lon=c14db$Longitude))
rownames(sites)  <- sites$id
sites  <- sites[,-1]
sp::coordinates(sites)  <- c('lon','lat')
sp::proj4string(sites)  <- sp::CRS("+proj=longlat +datum=WGS84")

# Compute distance matrix
d <- sp::spDists(sites,sites,longlat=TRUE)
# Compute spatial weights
w <- spweights(d,h=100)

# Define Breaks
breaks  <- BCADtoBP(c(-1000,-750,500,250,1,250))
timeRange  <- BCADtoBP(c(-1000,250))

