# Load Library ----
library(here)
library(dplyr)

# x  <- read.csv(here('data','prefecture_region_match.csv'))
# y  <- read.csv(here('data','prefectures_translations.csv'))
# xy  <- left_join(x,y,by=c('Prefecture'='Translation')) |> subset(!Region.x%in%c('Hokkaido','Okinawa'))
# xy  <- select(xy,PrefectureEn=Prefecture,JpNames,Region=Region.x,Region.koyama78=Region.alt,Area=Area.x,Surface=Area.y,Yayoi.start=Yayoi.start)
# write.csv(xy,row.names=FALSE,file=here('data','prefecture_data.csv'))

# Read External Data and Define Constants ----
# Settlement size ratio for Yayoi, as discussed in Koyama 1978
C.yayoi.haji  <- 1/3
# Koyama and Sugito 1984 Population Estimate
koyama84  <- read.csv(here('data','koyama_popestimate_1984.csv'))
koyama84  <- koyama84[-10,] #remove last row with total

# Number of Archaeological Sites from https://www.bunka.go.jp/seisaku/bunkazai/shokai/pdf/maizobunkazai.pdf, table 1, page 27
maizobunkazai2017  <- read.csv(here('data','maizobunkazai_2017.csv')) |> subset(!Prefecture%in%c('北海道','沖縄県','合計'))


# Lookup table for matching
pref.lookup  <- read.csv(here('data','prefecture_data.csv'))

# Prepare Data ----
maizobunkazai2017  <- left_join(maizobunkazai2017,pref.lookup,by=c('Prefecture'='JpNames'))
sitecounts  <- group_by(maizobunkazai2017,Region.koyama78) |> summarise(Yayoi.sites=sum(弥生),Haji.sites=sum(古墳)+sum(古代)) |> as.data.frame()
sitecounts  <- left_join(sitecounts,koyama84,by=c('Region.koyama78'='Region')) |> select(Region=Region.koyama78,Yayoi.pop84=Yayoi,Haji.pop=Haji,Yayoi.sites,Haji.sites)

sitecounts$PhSh  <- sitecounts$Haji.pop/sitecounts$Haji.sites

# Population Estimates for Koyama 1978 Regions ----

#Method 1: Koyama 1978
Ph.Sh.kanto  <- subset(sitecounts,Region=='Kanto')$PhSh
sitecounts$Yayoi.pop.est1  <- Ph.Sh.kanto * C.yayoi.haji * sitecounts$Yayoi.sites

#Method 2: Variable Sampling Fraction
sitecounts$Yayoi.pop.est2  <- sitecounts$PhSh * C.yayoi.haji * sitecounts$Yayoi.sites

#Method 3: Variable Sampling Fraction and Duration
prefect.counts  <- select(maizobunkazai2017,Prefecture=PrefectureEn,Region.koyama78,Yayoi.sites=弥生,Surface,Area,Yayoi.start)
prefect.counts  <- left_join(prefect.counts,select(sitecounts,Region,PhSh),by=c('Region.koyama78'='Region'))
prefect.counts$Yayoi.pop.est3  <- prefect.counts$PhSh * C.yayoi.haji * prefect.counts$Yayoi.sites * 1000/abs(prefect.counts$Yayoi.start - 250)
sitecounts  <- left_join(sitecounts,aggregate(Yayoi.pop.est3 ~ Region.koyama78, sum, data=prefect.counts),by=c('Region'='Region.koyama78'))

# barplot(t(as.matrix(sitecounts[,c(2,6,7,9)])),beside=T,names.arg=sitecounts$Region)

# Estimate based on Crema et al 2022 Regions ----
estimates  <- group_by(prefect.counts,Area) |> summarise(Surface=sum(Surface),EstPop=sum(Yayoi.pop.est3)) |> as.data.frame()
estimates$PopDens  <- estimates$EstPop/estimates$Surface





