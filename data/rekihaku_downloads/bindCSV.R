#### This scripts prepare the radiocarbon dates downloaded from the "Database of radiocarbon dates published in Japanese archaeological research reports) (URL: https://www.rekihaku.ac.jp/up-cgi/login.pl?p=param/esrd/db_param). 
#### Queries were carried out on the 4th of July 2023 by specifying "試料分類" (Material Classification) to "T: 陸産物" (Terrestial) and using C14 age intervals: 5000 - 3001, 3000 - 2001, 2000 - 1201, and 1200 - 100 to ensure the number of entries were below the download limit of 10,0000.  

library(here)
library(dplyr)
library(readxl)

c14.raw.1 <- read_excel(here('data','rekihaku_downloads','5000_3001_T.xlsx'))
c14.raw.2 <- read_excel(here('data','rekihaku_downloads','3000_2001_T.xlsx'))
c14.raw.3 <- read_excel(here('data','rekihaku_downloads','2000_1201_T.xlsx'))
c14.raw.4 <- read_excel(here('data','rekihaku_downloads','1200_100_T.xlsx'))

c14raw  <- rbind.data.frame(c14.raw.1,c14.raw.2,c14.raw.3,c14.raw.4) |> select(LabCode=試料番号,SamplingLocation=サンプル採取地点等) |> unique() |> subset(!is.na(LabCode) & !LabCode=='')
c14raw$LabCode  <- trimws(c14raw$LabCode)

c14raw <- aggregate(SamplingLocation~LabCode,toString,data=c14raw)

write.csv(c14raw,here('data','rekihaku_downloads','binded.csv'),row.names=FALSE)
