library(here)
library(dplyr)
load(here('data','sitedata.RData'))
load(here('data','c14data.RData'))

# Rice Arrival Dates ----
table1 = data.frame(regions=1:8,
			 m = c(-1039,-570,-910,-824,-648,-271,-152,-428),
			 hi = c(-1251,-735,-1061,-946,-754,-471,-434,-709),
			 lo = c(-872,-430,-779,-703,-560,-124,42,-203))


# Fujio settlement data ----
table1$site_fujio  <- table(sitedb$Area)

# Subset based on probability mass
foo  <- function(x1,y1,x2,y2)
{
	d  <- x1 - y1
	p  <- 1/d
	if (x1<=x2 & y1>=y2) {return(1)}
	if (x1>x2 & y1>x2) {return(0)}
	if (x1<y2 & y1<y2) {return(0)}
	if (x1>=x2 & y1>y2) {return(abs(x2-y1)*p)}
	if (x1<x2 & y1<y2) {return(abs(x1-y2)*p)}
	if (x1>=x2 & y1<=y2) {return(abs(x2-y2)*p)}
}

# Compute probabilities within window
sitedb$prob.within  <- unlist(apply(cbind(sitedb$start,sitedb$end,2900,1700),1,function(x){foo(x[1],x[2],x[3],x[4])}))
sitedb <- subset(sitedb,prob.within>0.5)
table1$site_fujio2  <- table(sitedb$Area)

# 14C data ----
table1$yayoidates  <- table(subset(c14db,Yayoi==TRUE)$RiceRegion) 
table1$d500dates  <- table(subset(c14db,D500==TRUE)$RiceRegion) 
table1$d750dates  <- table(subset(c14db,D750==TRUE)$RiceRegion) 
