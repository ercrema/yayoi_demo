library(here)
library(rcarbon)

# General Settings and Parameters:
anthropicGrep = c("住居","埋葬","竪穴建物","掘立柱","墓","包含層","土坑","ピット","土器","捨場","遺構","炉","人骨","木舟","住","柱","Pit","焼土","カマド","床面","溝中","溝底部","建物跡","木製品","埋土","水田","竪坑","羨道","集石","漆器")


# Read 14C Data ----
load(here('data','japanc14db_v03.3(210218).RData')) #version 3.0

# Cleaning ----
# Consider only dates with Retain=TRUE
c14db = subset(c14db,retain==TRUE)
# Add Booleans fields
c14db$anthropic = grepl(paste(anthropicGrep,collapse="|"),c14db$SamplingLocation) #sum(c14db$anthropic) 10683
save(c14db,file=here('data','c14data.RData'))
