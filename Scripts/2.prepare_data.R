library(readxl)
library("janitor")

# LOAD DATA -----

RR100s <- read.csv("Data/RR100s.csv")
RR100s <- RR100s[,-1]
RR100s$incub_date <- as.Date(as.character(RR100s$incub_date),"%y%m%d")
bdgp <- names(RR100s)

# CBA data
cba100lakes <- read_xlsx("Data/CBA_100Lakes_Master.xlsx")

# MERGE DATAFRAME -----

rr100lakes <- merge(cba100lakes, RR100s, by.y = "Sample_ID",by.x = "Lake_ID")
rr100lakes$lag_bdg <- difftime(rr100lakes$incub_date,rr100lakes$CBA_date,units = "days") %>% as.numeric()

rr100lakes$H <- 10^(-rr100lakes$pH)

rr100lakes$CN <- (rr100lakes$DOC/12.011)/(rr100lakes$DN/14.007)
rr100lakes$CP <- (rr100lakes$DOC/12.011)/(rr100lakes$DP*10^(-3)/30.97)
rr100lakes$NP <- (rr100lakes$DN/14.007)/(rr100lakes$DP*10^(-3)/30.97)
rr100lakes$CNP <-  (rr100lakes$DOC/12.011)/(rr100lakes$DN/14.007)/(rr100lakes$DP/30.97*10^6)

rr100lakes$s_350_400[which(rr100lakes$s_350_400 == 0)] <- NA
rr100lakes$SR <- rr100lakes$s_275_295/rr100lakes$s_350_400
rr100lakes$sVISa <- rr100lakes$a400/rr100lakes$DOC

#Altitude Langtjern
rr100lakes[grep("13763",rr100lakes$Sample_ID),which(names(rr100lakes) == "Altitude")] <- 516

#change names
names(rr100lakes)[which(names(rr100lakes)=="dmax")] <- "RR" # respiration rate, max speed of oxygen consumption
names(rr100lakes)[which(names(rr100lakes)=="hauc")] <- "OD" # oxygen demand = amount of O2 consumed
names(rr100lakes)[which(names(rr100lakes)=="hwidth")] <- "BdgT" # length of the growth phase
names(rr100lakes)[which(names(rr100lakes)=="TOT-P")] <- "TP" # length of the growth phase
names(rr100lakes)[which(names(rr100lakes)=="EC")] <- "EC_bio" # length of the growth phase
names(rr100lakes)[which(names(rr100lakes)=="pH")] <- "pH_bio" # length of the growth phase
names(rr100lakes)[which(names(rr100lakes)=="EC_Kje")] <- "EC" # length of the growth phase
names(rr100lakes)[which(names(rr100lakes)=="pH_Kje")] <- "pH" # length of the growth phase


write.csv(rr100lakes,"Data/rr100lakes.csv")
