library(here)
library(dplyr)
# Region Names ----
table1  <- data.frame(area=as.roman(1:8))

# Arrival Dates from Crema et al 2022, Model II https://doi.org/10.1126/sciadv.adc9171 ----
arrival <- read.csv('https://raw.githubusercontent.com/ercrema/yayoi_rice_dispersal/main/manuscript/supplementary_tables/table_S4.csv') |> subset(Model=='Model b')
table1$arrival.date  <- paste0(arrival[,3],' (',arrival[,4],'~',arrival[,5],')')

# 14C dates ----
load(here('data','c14data.RData'))
table1$d500dates  <- table(subset(c14db,D500==TRUE)$RiceRegion) 
table1$d750dates  <- table(subset(c14db,D750==TRUE)$RiceRegion) 
table1$yayoidates  <- table(subset(c14db,Yayoi==TRUE)$RiceRegion) 

# Fujio settlement data ----
load(here('data','sitedata.RData'))
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
table1$site_fujio  <- table(sitedb$Area)


# Population Estimate ----
pop.est  <- read.csv(here('results','pop_estimate_region.csv'))
table1$popdens <- pop.est$PopDens

# Define Column Names and Export ----
colnames(table1)  <- c('Region','Est. Arrival Date','N. of Dates (500yrs)','N. of Dates (750yrs)','N. of Dates (Yayoi period)','N. of sites','Pop. density est. (per km2')
write.csv(table1,here('tables','table1.csv'))
