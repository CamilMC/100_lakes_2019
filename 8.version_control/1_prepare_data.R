library(readxl)
library("janitor")

# LOAD DATA -----

RR100s <- read.csv("4.Respiration_rate/RR100s.csv")
RR100s <- RR100s[,-1]
RR100s$incub_date <- as.Date(as.character(RR100s$incub_date),"%y%m%d")
#RR100s$Sample_ID <- as.character(RR100s$Sample_ID)
bdgp <- names(RR100s)

# CBA data
cba100lakes <- read_xlsx("CBA100lakes_Master.xlsx")

#niva data
niva100 <- read.csv("7. 1000 lakes/niva1000.csv") %>% select(-c("dmax","hauc","hwidth","tmax","lag_bdg","Latitude","Longitude")) %>% distinct()
names(niva100)[which(names(niva100)=="Ca")] <- "Ca_NIVA"
names(niva100)[which(names(niva100)=="Na")] <- "Na_NIVA"
names(niva100)[which(names(niva100)=="K")] <- "K_NIVA"
names(niva100)[which(names(niva100)=="Mg")] <- "Mg_NIVA"
names(niva100)[which(names(niva100)=="DOC")] <- "DOC_NIVA"
names(niva100)[which(names(niva100)=="TOC.x")] <- "TOC_NIVA"
niva100$TOC.y <- NULL
niva100$sUVa <- NULL
niva100$incubation_date
#names(niva100)[which(names(niva100)=="Mg")] <- "Mg_NIVA"

missinglakes <- setdiff(lake$lake_id,niva100$Sample_ID) %>% as.data.frame()
emptylakes <-  missinglakes %>% cbind(matrix(data = NA, nrow = dim(missinglakes)[1], ncol = dim(niva100)[2]-1)) %>% setNames(names(niva100))

niva100 <- rbind(niva100, emptylakes)


# MERGE DATAFRAME -----

rr100lakes <- merge(cba100lakes, RR100s, by.y = "Sample_ID",by.x = "Lake_ID")
rr100lakes$lag_bdg <- difftime(rr100lakes$incub_date,rr100lakes$CBA_sample_date,units = "days") %>% as.numeric()

rr100lakes$DOC_umol <- rr100lakes$DOC/12 * 10^(3)
rr100lakes$hauc.DOCumol <- rr100lakes$hauc/rr100lakes$DOC_umol
rr100lakes$dmax.DOCumol <- rr100lakes$dmax/rr100lakes$DOC_umol

rr100lakes$dmax.hauc <- rr100lakes$dmax/rr100lakes$hauc

rr100lakes$H <- 10^(-rr100lakes$pH)

rr100lakes$CN <- rr100lakes$DOC_umol/(rr100lakes$DN/14.007*10^6)
rr100lakes$CP <- rr100lakes$DOC_umol/(rr100lakes$DP/30.97*10^6)

rr100lakes$s_350_400[which(rr100lakes$s_350_400 == 0)] <- 0.000001
rr100lakes$SR <- rr100lakes$s_275_295/rr100lakes$s_350_400

#Altitude Langtjern
rr100lakes[grep("13763",rr100lakes$Sample_ID),which(names(rr100lakes) == "Altitude")] <- 516

#change names
names(rr100lakes)[which(names(rr100lakes)=="dmax")] <- "Biodeg"
names(rr100lakes)[which(names(rr100lakes)=="hauc")] <- "BDOC"
names(rr100lakes)[which(names(rr100lakes)=="hwidth")] <- "BdgT"
names(rr100lakes)[which(names(rr100lakes)=="dmax.hauc")] <- "Biodeg.BDOC"
names(rr100lakes)[which(names(rr100lakes)=="dmax.DOCumol")] <- "Biodeg.DOC"
names(rr100lakes)[which(names(rr100lakes)=="LOC.DOC")] <- "BDOC.DOC"

write.csv(rr100lakes,"8.version_control/rr100lakes.csv")
write_xlsx(rr100lakes,"8.version_control/rr100lakes.xlsx")


# outliers -----
propTOC <- (all.pred$TOC-all.pred$abs254)
ggplot(all.pred,aes(x=abs254,label = Sample_ID))+geom_text(aes(y=DOC,col="DOC"))+geom_text(aes(y=TOC,col="TOC"))
ggplot(all.pred,aes(x=lake_name,y=propTOC))+geom_text(aes(label=Sample_ID))+geom_abline(slope=0,intercept= mean(propTOC),col="red")+
  geom_abline(slope=0,intercept= mean(propTOC)+sd(propTOC),col="pink")+
  geom_abline(slope=0,intercept= mean(propTOC)-sd(propTOC),col="pink")

# remove outliers ---
all.70 <- filter(all.pred,!Sample_ID %in% c("12777","13156","13196"))

# select and rename parameters ----

sel <- c("Sample_ID","hauc","tmax","hwidth","dmax","dmax.hauc","dmax.DOC","hauc.DOC","lag_bdg",
         "temp","long","lat","pH_lab","H","cond_µS.m","alk_meq.L",
         "TOC","TN","TP","DOC","DN","DP","CN","CP",
         "O2_µM", "CO2_µM","CH4_nM",
         "cells_counts.mL",
         "Fluoride_mg.L","Chloride_mg.L","Nitrite_mg.L","Bromide_mg.L","Sulfate_mg.L","Nitrate_mg.L",
         "Cr52","Ni58","Cu63","Zn64","As75","Cd114","Pb208","V51","Mn55","Co59",
         "abs254","abs410","SUVA","SAR","s_275_295","s_350_400","SR",
         "Altitude","Lake_Area","Ca","Mg","Na","K","Fe","Al","Al_R","Ca_NIVA","Na_NIVA")

lakes70 <- all.70[sel]
names(lakes70) <- c("Sample_ID","LOC","tmax","BdgT","Vmax","Vmax.LOC","Vmax.DOC","LOC.DOC","lag_bdg",
                    "temperature","long","lat","pH","H","conductivity","alkalinity",
                    "TOC","TN","TP","DOC","DN","DP","CN","CP",
                    "O2","CO2","CH4",
                    "cells",
                    "F","Cl","NO2","Br","SO4","NO3","Cr","Ni","Cu","Zn","As","Cd","Pb","V","Mn","Co",
                    "abs254","abs5410","SUVA","SAR","s_275_295","s_350_400","SR",
                    "Altitude","Lake_area","Ca","Mg","Na","K","Fe","Al","Al_NIVA","Ca_NIVA","Na_NIVA")
write.csv(lakes70,"5. 100 lakes/lakes70.csv")
