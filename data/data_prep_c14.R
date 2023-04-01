library(here)
library(rcarbon)
library(rnaturalearth)
library(sf)
library(maptools)
library(rgeos)
library(spdep)
source(here('src','dbscanID.R'))

# Read 14C Data ----
c14db  <- readRDS(here('data','c14db_0.2.1.Rds'))

# Subset to Anthropogeneic Dates ----
c14db <- subset(c14db,!is.na(Latitude)&!is.na(Longitude)&!PrefectureNameEn%in%c('Hokkaido','Okinawa'))
# NOTE : Any script for removing non-anthropic dates should be included here
# # Consider only dates with Retain=TRUE
# c14db = subset(c14db,retain==TRUE & !Region %in% c('Hokkaido','Okinawa'))
# # Add Booleans fields
# anthropicGrep = c("住居","埋葬","竪穴建物","掘立柱","墓","包含層","土坑","ピット","土器","捨場","遺構","炉","人骨","木舟","住","柱","Pit","焼土","カマド","床面","溝中","溝底部","建物跡","木製品","埋土","水田","竪坑","羨道","集石","漆器")
# c14db$anthropic = grepl(paste(anthropicGrep,collapse="|"),c14db$SamplingLocation) #sum(c14db$anthropic) 10683

# Aggregate into aritifical Sites using DBSCAN ----
# Add SiteID based on DBSCAN
c14db$SiteID  <- dbscanID(sitename=c14db$SiteNameEn,longitude = c14db$Longitude,latitude = c14db$Latitude,eps=100)

# Add Region Field and Predicted Arrival of Rice in each:
c14db$RiceRegion  <- NA
c14db$RiceRegion[which(c14db$PrefectureNameEn%in%c('Fukuoka','Saga','Nagasaki'))]  <- "I"
c14db$RiceRegion[which(c14db$PrefectureNameEn%in%c('Oita','Miyazaki','Kagoshima','Kumamoto'))]  <- "II"
c14db$RiceRegion[which(c14db$PrefectureNameEn%in%c('Ehime','Kochi')|c14db$Region=='Chugoku')]  <- "III"
c14db$RiceRegion[which(c14db$PrefectureNameEn%in%c('Kagawa','Tokushima')|c14db$Region=='Kansai')]  <- 'IV'
c14db$RiceRegion[which(c14db$PrefectureNameEn%in%c('Kanagawa')|c14db$Region=='Chubu')]  <- "V"
c14db$RiceRegion[which(c14db$PrefectureNameEn%in%c('Saitama','Gunma','Tokyo','Chiba','Tochigi','Ibaraki'))]  <- "VI"
c14db$RiceRegion[which(c14db$PrefectureNameEn%in%c('Yamagata','Fukushima','Miyagi','Iwate','Akita'))]  <- "VII"
c14db$RiceRegion[which(c14db$PrefectureNameEn%in%c('Aomori'))]  <- "VIII"
# any(is.na(c14db$RiceRegion))
# c14db$PrefectureNameEn[which(is.na(c14db$RiceRegion))]

# Add rice arrival date in BP
c14db$ricearrival <- NA
c14db$ricearrival[c14db$RiceRegion=='I']  <- BCADtoBP(-1039)
c14db$ricearrival[c14db$RiceRegion=='II']  <- BCADtoBP(-570)
c14db$ricearrival[c14db$RiceRegion=='III']  <- BCADtoBP(-910)
c14db$ricearrival[c14db$RiceRegion=='IV']  <- BCADtoBP(-824)
c14db$ricearrival[c14db$RiceRegion=='V']  <- BCADtoBP(-648)
c14db$ricearrival[c14db$RiceRegion=='VI']  <- BCADtoBP(-271)
c14db$ricearrival[c14db$RiceRegion=='VII']  <- BCADtoBP(-152)
c14db$ricearrival[c14db$RiceRegion=='VIII']  <- BCADtoBP(-428)

# Handle Dates
c14db$C14Age = c14db$UnroundedCRA
i = which(is.na(c14db$C14Age))
c14db$C14Age[i] = c14db$CRA[i]

c14db$C14Error = c14db$UnroundedCRAError
i = which(is.na(c14db$C14Error))
c14db$C14Error[i] = c14db$CRAError[i]

c14db  <- subset(c14db,!is.na(C14Age) & !is.na(C14Error))

# Spatial Data -----
# Prefecture Based
win <- ne_states(country = "japan",returnclass = 'sf') |> subset(!name_vi %in% c("Okinawa","Hokkaido"))
Npref <- nrow(win)
win$ID  <-  1:Npref
W_nb <- poly2nb(win, row.names =  win$ID)
nbInfo <- nb2WB(W_nb)

# RiceRegion Based
win$riceregion <- NA
win$riceregion[which(win$name%in%c('Fukuoka','Saga','Nagasaki'))]  <- "I"
win$riceregion[which(win$name%in%c('Oita','Miyazaki','Kagoshima','Kumamoto'))]  <- "II"
win$riceregion[which(win$name%in%c('Ehime','Kochi')|win$region=='Chugoku')]  <- "III"
win$riceregion[which(win$name%in%c('Kagawa','Tokushima')|win$region=='Kinki')]  <- 'IV'
win$riceregion[which(win$name%in%c('Kanagawa')|win$region=='Chubu')]  <- "V"
win$riceregion[which(win$name%in%c('Saitama','Gunma','Tokyo','Chiba','Tochigi','Ibaraki'))]  <- "VI"
win$riceregion[which(win$name%in%c('Yamagata','Fukushima','Miyagi','Iwate','Akita'))]  <- "VII"
win$riceregion[which(win$name%in%c('Aomori'))]  <- "VIII"
win.riceregion  <- group_by(win,riceregion) |> summarise()

Nregions  <- nrow(win.riceregion)
win.riceregion$ID  <- 1:Nregions
W_nb.rice <- poly2nb(win.riceregion, row.names =  win.riceregion$ID)
nbInfo.rice <- nb2WB(W_nb.rice)

# Store everything in a R image file ----

save(c14db,win,win.riceregion,Npref,Nregions,W_nb,nbInfo,W_nb.rice,nbInfo.rice,file=here('data','c14data.RData'))
