library(here)
library(rcarbon)

# Read 14C Data ----
load(here('data','japanc14db_v04.0(210715).RData')) #version 3.0

# Preprocessing ----
# Consider only dates with Retain=TRUE
c14db = subset(c14db,retain==TRUE & !Region %in% c('Hokkaido','Okinawa'))
# Add Booleans fields
anthropicGrep = c("住居","埋葬","竪穴建物","掘立柱","墓","包含層","土坑","ピット","土器","捨場","遺構","炉","人骨","木舟","住","柱","Pit","焼土","カマド","床面","溝中","溝底部","建物跡","木製品","埋土","水田","竪坑","羨道","集石","漆器")
c14db$anthropic = grepl(paste(anthropicGrep,collapse="|"),c14db$SamplingLocation) #sum(c14db$anthropic) 10683
# Add SiteID based on DBSCAN
source(here('R','dbscanID.R'))
c14db$SiteID=dbscanID(sitename=c14db$SiteName,longitude = c14db$Longitude,latitude = c14db$Latitude,eps=0.1)
# Add Region ID
c14db$RegionID <- as.numeric(factor(c14db$Region, levels = c("Kyushu", "Chugoku","Shikoku", "Kansai", "Chubu", "Kanto", "Tohoku")))
c14db$Region2 <- c14db$Region
c14db$Region2[c14db$Region2 %in% c("Shikoku","Chugoku")] = "Chugoku-Shikoku"
c14db$RegionID2 <- as.numeric(factor(c14db$Region2, levels = c("Kyushu", "Chugoku-Shikoku", "Kansai", "Chubu", "Kanto", "Tohoku"),
ordered = T))
c14db$Region3 <- c14db$Region2
c14db$Region3[c14db$PrefectureNameEn %in% c("Aomori","Akita","Iwate")] = "NorthernTohoku"
c14db$Region3[c14db$PrefectureNameEn %in% c("Yamagata","Miyagi","Fukushima","Niigata")] = "SouthernTohoku"
c14db$RegionID3 <- as.numeric(factor(c14db$Region3, levels = c("Kyushu","Chugoku-Shikoku","Kansai","Chubu","Kanto","SouthernTohoku","NorthernTohoku"),ordered = T))


c14db$C14Age = c14db$UnroundedCRA
i = which(is.na(c14db$C14Age))
c14db$C14Age[i] = c14db$CRA[i]

c14db$C14Error = c14db$UnroundedCRAError
i = which(is.na(c14db$C14Error))
c14db$C14Error[i] = c14db$CRAError[i]




save(c14db,file=here('data','c14data.RData'))
