library(readxl)
library("janitor")
# LOAD DATA -----

RR100s <- read.csv("4. Respiration rate/RR100s.csv")
RR100s <- RR100s[,-1]
RR100s$incub_date <- as.Date(as.character(RR100s$incub_date),"%y%m%d")
#RR100s$Sample_ID <- as.character(RR100s$Sample_ID)
bdgp <- names(RR100s)

#TOC data
om_data <- read_xlsx("5. 100 lakes/om_data.xlsx")

#sampling data
lake_data <- read_xlsx("5. 100 lakes/100LakesData.xlsx",col_types = c("text","text","numeric","numeric","numeric","numeric","numeric","date","numeric"))
lake <- lake_data[!duplicated(lake_data$lake_id),]

#absorbance data
abs_data <- read_xlsx("5. 100 lakes/100lakes_all_absorbances.xlsx")
abs <- abs_data %>% filter(WL.nm==254) %>% t() %>% as.data.frame() %>% tibble::rownames_to_column() %>% filter(rowname != "WL.nm") %>% setNames(c("lake_ID","abs254"))
abs$abs254 <- as.numeric(abs$abs254)
abs$abs410 <- abs_data %>% filter(WL.nm==410) %>% t() %>% row_to_names(row_number = 1) %>% as.numeric()

# spectral slopes

s_275_295 <- getslope(abs_data,275,295)
s_350_400 <- getslope(abs_data,350,400)

#bacterial data

dna <- read_xlsx("5. 100 lakes/100lakes_DNA.xlsx", col_types = c("text","text","numeric","text","numeric","numeric","numeric"))

#gas data
conc_gas <- read_xlsx("5. 100 lakes/100Lakes_gases_conc_rates.xlsx",col_types = c("text","skip","skip","skip","skip","skip","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric","numeric"))

#water chemistry data
chem_data <- read_xlsx("5. 100 lakes/100lakes_chemistry.xlsx",col_types = c("skip","skip","text","skip","numeric","numeric","numeric","numeric","numeric","skip"))

#anion data
anion_data <- read_xlsx("5. 100 lakes/100lakes_anions.xlsx",col_types = c("skip","text","numeric","numeric","numeric","numeric","numeric","numeric","numeric")) 

# trace metals
tm_data <- read_xlsx("5. 100 lakes/100lakes_ICPMS.xlsx",sheet = "Sheet3")

#cations
bcations <- read_xlsx("5. 100 lakes/100lakes_base_cations.xlsx",col_types = c("text","numeric","numeric","numeric","numeric","numeric","numeric"))

#niva data
niva100 <- read.csv("7. 1000 lakes/niva1000.csv") %>% 
  select(c("Sample_ID","Lake_Area","Altitude","Ca","Mg","Na","K","Al_R","Al_Il","NH4.N")) %>% 
  distinct()
names(niva100) <- c("Sample_ID","Lake_Area","Altitude","Ca_NIVA","Mg_NIVA","Na_NIVA","K_NIVA","Al_R","Al_Il","NH4.N")
missinglakes <- setdiff(lake$lake_id,niva100$Sample_ID) %>% as.data.frame()
emptylakes <-  missinglakes %>% cbind(matrix(data = NA, nrow = dim(missinglakes)[1], ncol = dim(niva100)[2]-1)) %>% setNames(names(niva100))

niva100 <- rbind(niva100, emptylakes)

#MPAES
mpaes <- read_xlsx("5. 100 lakes/100lakes_MPAES.xlsx",sheet = "Summary_MPAES",
                   col_types = "numeric") %>%
      select(c("Sample","Ca","Na","Mg","K","Al","Fe"))


# MERGE DATAFRAME - all.pred -> all.parameters.csv -----

all.pred <- RR100s[,bdgp] %>% merge(lake,by.x="Sample_ID",by.y="lake_id") %>%
  merge(chem_data,by.x="Sample_ID",by.y="Lake_ID") %>% 
  merge(om_data,by.x="Sample_ID",by.y="ID") %>% 
  merge(abs,by.x="Sample_ID",by.y="lake_ID")%>%
  merge(s_275_295,by.x = "Sample_ID", by.y = "Lake_ID") %>% 
  merge(s_350_400, by.x = "Sample_ID", by.y = "Lake_ID") %>% 
  merge(dna, by.x = "Sample_ID", by.y = "Lake_ID") %>%
  merge(conc_gas, by.x = "Sample_ID", by.y = "Lake_ID") %>% merge(anion_data,by.x = "Sample_ID",by.y="ID") %>% 
  merge(tm_data,by.x = "Sample_ID",by.y = "ID") %>%
  merge(bcations,by.x = "Sample_ID", by.y = "Sample_ID") %>% 
  merge(niva100, by.x = "Sample_ID", by.y = "Sample_ID") %>%
  merge(mpaes, by.x = "Sample_ID", by.y = "Sample")

all.pred$lag_bdg <- difftime(all.pred$incub_date,all.pred$date,units = "days") %>% as.numeric()
all.pred$SUVA <- all.pred$abs254/(all.pred$TOC*0.90)
all.pred$abs410[all.pred$abs410 <= 0 ] <- NA
all.pred$SAR <- all.pred$abs254/all.pred$abs410
all.pred$s_350_400[all.pred$s_350_400 <= 0 ] <- NA
all.pred$SR <- all.pred$s_275_295/all.pred$s_350_400

all.pred$hauc.TOC <- all.pred$hauc/all.pred$TOC
all.pred$auc.TOC <- all.pred$auc/all.pred$TOC
all.pred$dmax.TOC <- all.pred$dmax/all.pred$TOC
all.pred$tmax.TOC <- all.pred$tmax / all.pred$TOC
all.pred$hwidth.TOC <- all.pred$hwidth/ all.pred$TOC

all.pred$hauc.DOC <- all.pred$hauc/all.pred$DOC
all.pred$auc.DOC <- all.pred$auc/all.pred$DOC
all.pred$dmax.DOC <- all.pred$dmax/all.pred$DOC
all.pred$tmax.DOC <- all.pred$tmax / all.pred$DOC
all.pred$hwidth.DOC <- all.pred$hwidth/ all.pred$DOC

all.pred$DOC_umol.L <- all.pred$DOC/12 * 10^(3)
all.pred$hauc.DOCumol <- all.pred$hauc/all.pred$DOC_umol.L
all.pred$auc.DOCumol <- all.pred$auc/all.pred$DOC_umol.L
all.pred$dmax.DOCumol <- all.pred$dmax/all.pred$DOC_umol.L
all.pred$tmax.DOCumol <- all.pred$tmax / all.pred$DOC_umol.L
all.pred$hwidth.DOCumol <- all.pred$hwidth/ all.pred$DOC_umol.L

all.pred$dmax.hauc <- all.pred$dmax/all.pred$hauc

all.pred$H <- 10^(-all.pred$pH_lab)

all.pred$CN <- all.pred$DOC/all.pred$DN
all.pred$CP <- all.pred$DOC/all.pred$DP

#Altitude Langtjern
all.pred[grep("13763",all.pred$Sample_ID),which(names(all.pred) == "Altitude")] <- 516

write.csv(all.pred,"5. 100 lakes/all.parameters.csv")

# SELECT DATA: waterchem -----

all.pred <- read.csv("5. 100 lakes/all.parameters.csv")

#removes duplicates (see "200330_duplicates") -----
sel <- c("Sample_ID","hauc","tmax","hwidth","dmax","dmax.hauc","dmax.DOCumol","hauc.DOCumol","lag_bdg", "temp","long","lat","week","pH_lab","H",
         "cond_µS.m","alk_meq.L",
         "TOC","TN","TP","DOC","DN","DP",
         "O2_µM", "N2_µM","CO2_µM","CH4_nM","N2O_nM","O2_nM.h","N2_nM.h","CO2_nM.h","CH4_nM.h","N2O_nM.h",
         "NVE_number","DNA_ug.mL","filtered_mL","cells_counts.mL",
         "Fluoride_mg.L","Chloride_mg.L","Nitrite_mg.L","Bromide_mg.L","Sulfate_mg.L","Nitrate_mg.L",
         "Cr52","Ni58","Cu63","Zn64","As75","Cd114","Pb208","V51","Mn55","Co59",
         "abs254","abs410","SUVA","SAR","s_275_295","s_350_400","SR",
         "Altitude","Lake_Area","Ca","Mg","Na","K","Fe","Al","Al_R","Al_Il","Ca_NIVA","Na_NIVA")
waterchem <- all.pred[sel]
waterchem$Sample_ID <- as.numeric(waterchem$Sample_ID)
write.csv(waterchem,"5. 100 lakes/waterchem.csv")

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

