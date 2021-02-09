#libraries -----

library("sf")
library("ggrepel")
library("dplyr")
library(writexl)
library(readxl)
library(ggplot2)

library(dendextend)
library("effectsize")
library(factoextra)
library(viridisLite)
library(reshape2)

library(mice)
library(glmnet)

source("8.version_control/100lakes_analysis_functions.R")
norgemap <- st_read("3.Maps/TM_WORLD_BORDERS-0.3.shp")
source("8.version_control/bacteria.R")

# -----
# load the world map -----

library("rnaturalearth")
library(rgeos)

world <- ne_countries(scale = "medium", returnclass = "sf")

# -----
# load data -----

lakes73 <- read.csv("8.version_control/rr100lakes.csv")
lakes73$SARvis[which(lakes73$SARvis == "Inf")] <- "NA"
# ----
# summary stats RR, BdgT and OD -----

cor.test(lakes73$RR,lakes73$OD)[c("estimate","p.value")] 
cor.test(lakes73$BdgT,lakes73$OD)[c("estimate","p.value")] 
cor.test(lakes73$RR,lakes73$BdgT)[c("estimate","p.value")] 


# ----
# separation of eutrophic/oligotrophic/dystrophic lakes ------
lakes73$DNOM <- NA
for(i in 1:length(lakes73$TOC)){
  if (lakes73$TOC[i] >= mean(lakes73$TOC)){
    lakes73$DNOM[i] <- 1
    } else {
    lakes73$DNOM[i] <- 0
    }
}

lakes73$nuts <- NA
for(i in 1:length(lakes73$NP)){
  if ((lakes73$NP[i]) > (mean(lakes73$NP))){
    lakes73$nuts[i] <- 1
  } else {
    lakes73$nuts[i] <- 0
  }
}

lakes73$nuts <- NA
for(i in 1:length(lakes73$NP)){
  if ((lakes73$TN[i]/lakes73$TP[i]) >= (mean(lakes73$TN/lakes73$TP))){
    lakes73$nuts[i] <- 1
  } else {
    lakes73$nuts[i] <- 0
  }
}

lakes73$trophic <- NA
for (i in 1:length(lakes73$Lake_ID)){
  if(lakes73$DNOM[i] == 0 & lakes73$nuts[i] == 0){
    lakes73$trophic[i] <- "eu"
  }else if(lakes73$DNOM[i] == 0 & lakes73$nuts[i] == 1){
  lakes73$trophic[i] <- "oligo"
  } else if(lakes73$DNOM[i] ==1 & lakes73$nuts[i] == 0){
  lakes73$trophic[i] <- "dys" 
  } else {
    lakes73$trophic[i] <- "eu"
  }
}

ggplot()+geom_sf(data=norgemap)+
  geom_point(data=lakes73,aes(x=Long,y=Lat,col=trophic),size=6)+
  theme_void(base_size = 26)+labs(x="",y="",col="trophic state")+
  scale_color_manual(values=colors_cluster)+
  theme(legend.position = c(.8,.15))+
  xlim(5,15)+ylim(58,62)

ggplot()+geom_sf(data=norgemap)+
  geom_point(data=filter(lakes73,RR > 10),aes(x=Long,y=Lat),col = "red",size = 5)+
  geom_point(data=lakes73,aes(x=Long,y=Lat,col=trophic),size=4)+
  theme_void(base_size = 26)+labs(x="",y="",col="trophic")+
  scale_color_manual(values=colors_cluster)+
  theme(legend.position = c(.75,.15))+
  xlim(5,15)+ylim(58,63)

# -----
# alternative sorting of trophic status -----

lakes73$CNP <- lakes73$TOC/(lakes73$TN/lakes73$TP)

ggplot(lakes73)+geom_boxplot(aes(x="CNP",y=TOC/(TN/TP)),outlier.shape = NA)+geom_text_repel(aes(x="CNP",y=TOC/(TN/TP),label=Lake_ID))


lakes73$trophic <- NA
for (i in 1:length(lakes73$Lake_ID)){
  if(lakes73$CNP[i] > quantile(lakes73$CNP)[[4]]){
    lakes73$trophic[i] <- "eu"
  }else if(lakes73$CNP[i] < quantile(lakes73$CNP)[[3]]){
    lakes73$trophic[i] <- "oligo"
  } else {
    lakes73$trophic[i] <- "dys" 
    }
}

ggplot()+geom_sf(data=norgemap)+
  geom_point(data=lakes73,aes(x=Long,y=Lat,col=trophic),size=6)+
  theme_void(base_size = 26)+labs(x="",y="",col="trophic state")+
  scale_color_manual(values=colors_cluster)+
  theme(legend.position = c(.8,.15))+
  xlim(5,15)+ylim(58,62)

# -----
# alternative sorting of trophic status by P-----

lakes73$DNOM <- NA
for(i in 1:length(lakes73$TOC)){
  if (lakes73$TOC[i] >= mean(lakes73$TOC)){
    lakes73$DNOM[i] <- 1
  } else {
    lakes73$DNOM[i] <- 0
  }
}

lakes73$nuts <- NA
for(i in 1:length(lakes73$TP)){
  if ((lakes73$TP[i]) > (mean(lakes73$TP))){
    lakes73$nuts[i] <- 1
  } else {
    lakes73$nuts[i] <- 0
  }
}

lakes73$trophic <- NA
for (i in 1:length(lakes73$Lake_ID)){
  if(lakes73$DNOM[i] == 0 & lakes73$nuts[i] == 1){
    lakes73$trophic[i] <- "eu"
  }else if(lakes73$DNOM[i] == 0 & lakes73$nuts[i] == 0){
    lakes73$trophic[i] <- "oligo"
  } else if(lakes73$DNOM[i] ==1 & lakes73$nuts[i] == 0){
    lakes73$trophic[i] <- "dys" 
  } else {
    lakes73$trophic[i] <- "eu"
  }
}

ggplot()+geom_sf(data=norgemap)+
  geom_point(data=lakes73,aes(x=Long,y=Lat,col=trophic),size=6)+
  theme_void(base_size = 26)+labs(x="",y="",col="trophic state")+
  scale_color_manual(values=c("sienna","chartreuse4","skyblue3"))+
  theme(legend.position = c(.8,.15))+
  xlim(5,15)+ylim(58,62)

# -----
# separation of eutrophic/oligotrophic/dystrophic lakes based on Nünberg (1996) -----
lakes73$DNOM <- NA
for(i in 1:length(lakes73$DOC)){
  if (lakes73$DOC[i] > mean(lakes73$DOC)){
    lakes73$DNOM[i] <- 1
  } else {
    lakes73$DNOM[i] <- 0
  }
}

lakes73$indexN <- NA

for(i in 1:length(lakes73$TN)){

    if (lakes73$TN[i] < 0.35){
    lakes73$indexN[i] <- "o" 
    } else if (lakes73$TN[i] >= 0.35 & lakes73$TN[i] < 0.65){
    lakes73$indexN[i] <- "m"
    } else if (lakes73$TN[i] >= 0.65 & lakes73$TN[i] < 1.2){
    lakes73$indexN[i] <- "e"
    } else if (lakes73$TN[i] >= 1.2){
    lakes73$indexN[i] <- "he"
    }
}

lakes73$indexP <- NA

for(i in 1:length(lakes73$TP)){
  
  if (lakes73$TP[i] < 10){
    lakes73$indexP[i] <- "o" 
  } else if (lakes73$TP[i] >= 10 & lakes73$TP[i] < 30){
    lakes73$indexP[i] <- "m"
  } else if (lakes73$TP[i] >= 30 & lakes73$TP[i] < 100){
    lakes73$indexP[i] <- "e"
  } else if (lakes73$TP[i] >= 100){
    lakes73$indexP[i] <- "he"
  }
}


lakes73$trophic <- NA
for (i in 1:length(lakes73$Lake_ID)){
  if(lakes73$indexN[i] == "o" & lakes73$indexP[i] == "o"){
    lakes73$trophic[i] <- "o"
  }else if(lakes73$indexN[i] == "m" & lakes73$indexP[i] == "m"){
    lakes73$trophic[i] <- "m"
  } else if(lakes73$indexN[i] == "e" & lakes73$indexP[i] == "e"){
    lakes73$trophic[i] <- "e" 
  } else if(lakes73$indexN[i] == "he" & lakes73$indexP[i] == "he"){
    lakes73$trophic[i] <- "he" 
  } else {
    lakes73$trophic[i] <- "tbd"
  }
}

ggplot()+geom_sf(data=norgemap)+
  geom_point(data=lakes73,aes(x=Long,y=Lat,col=indexN),size=6)+
  theme_void(base_size = 26)+labs(x="",y="",col="trophicP")+
  scale_color_manual(values=colors_cluster)+
  theme(legend.position = c(.75,.15))+
  xlim(5,15)+ylim(58,63)

# -----
# dataframes with less parameters ------

lakes73num <- lakes73 %>% sapply(as.numeric) %>% as.data.frame() %>% select_if(not_all_na)

lakes65 <- filter(lakes73,RR<10)
lakes65num <- lakes65 %>% sapply(as.numeric) %>% as.data.frame() %>% select_if(not_all_na)

lakes8 <- filter(lakes73,RR>10)
lakes8num <- lakes8 %>% sapply(as.numeric) %>% as.data.frame() %>% select_if(not_any_na)

# -----
# distribution of data ----
pdf("8.version_control/distribution.pdf")
hist.df(lakes73num)
dev.off()

# -----
# pca with all lakes -----

pdf("8.version_control/pca_plots.pdf",width=10,height = 8)

lakes73num$trophic <- NULL
lakes73pca <- lakes73num %>% fill.na() %>% select(c("RR","OD","BdgT","CN","TOC","CP","DOC","tmax","ox_initial","EC_Kje","DN","TN","NO2","NO3","K","Ca","pH_Kje","Alkalinity","c_O2","Mg"))
#%>% select(c("DN","DP","DOC","CN","CP","SUVA","SR","Ca","Fe","H","RR","OD","BdgT","O2","CO2","cells_counts.mL"))
rownames(lakes73pca) <- lakes73$Lake_ID
pca73lakes <- prcomp(~.,data=lakes73pca,scale.=T,center=T)
pcafviz <- fviz_pca_biplot(pca73lakes, repel=T,label = c("var"),col.ind = "contrib",
                           gradient.cols =  viridis(3,end=0.8,direction = -1),arrowsize=1,labelsize=6,title="")+ theme_light(base_size=24) 
pcafviz23 <- fviz_pca_biplot(pca73lakes, axes = c(2,3),repel=T,label = c("var"),col.ind = "contrib",
                             gradient.cols =  viridis(3,end=0.8,direction = -1),arrowsize=1,labelsize=6,title="")+ theme_light(base_size=24) 
plot(pcafviz)
plot(pcafviz23)

# -----
# pca without all outliers  -----

lakes65pca <- lakes65 %>% fill.na() %>% select(c("RR","OD","BdgT","lag_bdg","CN","K","Na","EC_Kje","Mn","Mg","As","c_O2","pH","tmax","DN","NO2","TN","Ca","Alkalinity","TP")) %>% sapply(as.numeric) %>% as.data.frame()

rownames(lakes65pca) <- lakes65$Lake_ID
pca65lakes <- prcomp(~.,data=lakes65pca,scale.=T,center=T)
pcafviz <- fviz_pca_biplot(pca65lakes, repel=T,label = c("var"),col.ind = "contrib",
                           gradient.cols =  viridis(3,end=0.8,direction = -1),arrowsize=1,labelsize=6,title="")+ theme_light(base_size=24) 
pcafviz23 <- fviz_pca_biplot(pca65lakes, axes = c(2,3),repel=T,label = c("var"),col.ind = "contrib",
                             gradient.cols =  viridis(3,end=0.8,direction = -1),arrowsize=1,labelsize=6,title="")+ theme_light(base_size=24) 
plot(pcafviz)
plot(pcafviz23)
# -----
# pca of outliers -----
lakes8pca <- lakes73num %>% fill.na() %>% filter(RR > 10) %>% select(c("RR","OD","BdgT","SARuv","Zn","tmax","TOC","DOC","CN","CP","V","Cr"))

pca8lakes <- prcomp(~.,data=lakes8pca,scale.=T,center=T)
pcafviz <- fviz_pca_biplot(pca8lakes, repel=T,label = "var",col.var = "contrib",col.ind = "contrib",
                           select.var = list(contrib = 15), 
                           gradient.cols =  viridis(3,end=0.8,direction = -1),arrowsize=1,labelsize=6,title="")+ theme_light(base_size=24) 
plot(pcafviz)

pcafviz23 <- fviz_pca_biplot(pca8lakes, axes = c(2,3),repel=T,label = "var",col.var = "contrib",col.ind = "contrib",
                             select.var = list(contrib = 15),
                             gradient.cols =  viridis(3,end=0.8,direction = -1),arrowsize=1,labelsize=6,title="")+ theme_light(base_size=24) 
plot(pcafviz23)

dev.off()
# -----
# clean PCA ----

pdf("8.version_control/pca_plots.pdf",width=10,height = 8)

lakes73pca <- lakes73num %>% fill.na() %>% select(c("RR","OD","BdgT","CN","TOC","CP","EC_Kje","TN","pH_Kje","c_O2"))
rownames(lakes73pca) <- lakes73$Lake_ID
pca73lakes <- prcomp(~.,data=lakes73pca,scale.=T,center=T)
pcafviz <- fviz_pca_biplot(pca73lakes, repel=T,label = c("var"),col.ind = "contrib",
                           gradient.cols =  viridis(3,end=0.8,direction = -1),arrowsize=1,labelsize=6,title="")+ 
          theme_light(base_size=24)+
          labs(title = "73 lakes")

pcafviz23 <- fviz_pca_biplot(pca73lakes, axes = c(2,3),repel=T,label = c("var"),col.ind = "contrib",
                             gradient.cols =  viridis(3,end=0.8,direction = -1),arrowsize=1,labelsize=6,title="")+
            theme_light(base_size=24)+
            labs(title = "73 lakes")

plot(pcafviz)
plot(pcafviz23)
# -----
# pca without all outliers  -----
lakes65 <- filter(lakes73,RR<10)
lakes65pca <- lakes65 %>% fill.na() %>% select(c("RR","OD","BdgT","CN","EC_Kje","c_O2","p_CO2","pH_Kje","TN","TP")) %>% sapply(as.numeric) %>% as.data.frame()

rownames(lakes65pca) <- lakes65$Lake_ID
pca65lakes <- prcomp(~.,data=lakes65pca,scale.=T,center=T)
pcafviz <- fviz_pca_biplot(pca65lakes, repel=T,label = c("var"),col.ind = "contrib",
                           gradient.cols =  viridis(3,end=0.8,direction = -1),arrowsize=1,labelsize=6,title="")+ 
          theme_light(base_size=24)+
          labs(title = "Without outliers")
pcafviz23 <- fviz_pca_biplot(pca65lakes, axes = c(2,3),repel=T,label = c("var"),col.ind = "contrib",
                             gradient.cols =  viridis(3,end=0.8,direction = -1),arrowsize=1,labelsize=6,title="")+ 
          theme_light(base_size=24)+
          labs(title = "Without outliers")

plot(pcafviz)
plot(pcafviz23)

# -----
# pca of outliers -----
lakes8pca <- lakes73num %>% fill.na() %>% filter(RR > 10) %>% select(c("RR","OD","BdgT","SARuv","Zn","TOC","CN","CP","V","Cr"))

pca8lakes <- prcomp(~.,data=lakes8pca,scale.=T,center=T)
pcafviz <- fviz_pca_biplot(pca8lakes, repel=T,label = "var",col.ind = "contrib",
                           select.var = list(contrib = 15), 
                           gradient.cols =  viridis(3,end=0.8,direction = -1),arrowsize=1,labelsize=6,title="")+ 
          theme_light(base_size=24)+
          labs(title = "outliers")
plot(pcafviz)

pcafviz23 <- fviz_pca_biplot(pca8lakes, axes = c(2,3),repel=T,label = "var",col.ind = "contrib",
                             select.var = list(contrib = 15),
                             gradient.cols =  viridis(3,end=0.8,direction = -1),arrowsize=1,labelsize=6,title="")+ 
            theme_light(base_size=24)+
            labs(title = "outliers")
plot(pcafviz23)

dev.off()


# -----
# pca of DNOM characteristics -----
dnom_ctq <- c("RR","BdgT","SUVA","SARuv","pH","DOC","TN","TP")

pdf("8.version_control/pca_dnom_characteristics.pdf",width=10,height = 8)

###
dnom_ctq73 <- lakes73num %>% select(all_of(dnom_ctq))
pca_dnom_ctq <- prcomp(~.,data=dnom_ctq73,scale.=T,center=T)
pcafviz <- fviz_pca_biplot(pca_dnom_ctq, repel=T,label = "var",col.ind = "contrib",
                           select.var = list(contrib = 15), 
                           gradient.cols =  viridis(3,end=0.8,direction = -1),arrowsize=1,labelsize=6,title="")+ 
  theme_light(base_size=24)+
  labs(title = "all lakes (73)")
plot(pcafviz)

pcafviz23 <- fviz_pca_biplot(pca_dnom_ctq, axes = c(2,3),repel=T,label = "var",col.ind = "contrib",
                             select.var = list(contrib = 15), 
                             gradient.cols =  viridis(3,end=0.8,direction = -1),arrowsize=1,labelsize=6,title="")+ 
  theme_light(base_size=24)+
  labs(title = "all lakes (73)")
plot(pcafviz23)

###
dnom_ctq65 <- lakes65num %>% select(all_of(dnom_ctq))
pca_dnom_ctq <- prcomp(~.,data=dnom_ctq65,scale.=T,center=T)
pcafviz <- fviz_pca_biplot(pca_dnom_ctq, repel=T,label = "var",col.ind = "contrib",
                           select.var = list(contrib = 15), 
                           gradient.cols =  viridis(3,end=0.8,direction = -1),arrowsize=1,labelsize=6,title="")+ 
  theme_light(base_size=24)+
  labs(title = "main lakes (65")
plot(pcafviz)

pcafviz23 <- fviz_pca_biplot(pca_dnom_ctq, axes = c(2,3),repel=T,label = "var",col.ind = "contrib",
                             select.var = list(contrib = 15), 
                             gradient.cols =  viridis(3,end=0.8,direction = -1),arrowsize=1,labelsize=6,title="")+ 
  theme_light(base_size=24)+
  labs(title = "main lakes (65)")
plot(pcafviz23)

###
dnom_ctq8 <- lakes8num %>% select(all_of(dnom_ctq))
pca_dnom_ctq <- prcomp(~.,data=dnom_ctq8,scale.=T,center=T)
pcafviz <- fviz_pca_biplot(pca_dnom_ctq, repel=T,label = "var",col.ind = "contrib",
                           select.var = list(contrib = 15), 
                           gradient.cols =  viridis(3,end=0.8,direction = -1),arrowsize=1,labelsize=6,title="")+ 
  theme_light(base_size=24)+
  labs(title = "outliers (8)")
plot(pcafviz)

pcafviz23 <- fviz_pca_biplot(pca_dnom_ctq, axes = c(2,3),repel=T,label = "var",col.ind = "contrib",
                           select.var = list(contrib = 15), 
                           gradient.cols =  viridis(3,end=0.8,direction = -1),arrowsize=1,labelsize=6,title="")+ 
  theme_light(base_size=24)+
  labs(title = "outliers (8)")
plot(pcafviz23)


lakes73corpca <- lakes73cor %>% select(c("RR.DOC","SUVA","SARuv","pH","DOC","TN","TP"))

pcalakes73cor <- prcomp(~., data = lakes73corpca,scale.=T,center=T)
pcafviz <- fviz_pca_biplot(pcalakes73cor, repel=T,label = "var",col.ind = "contrib",
                           select.var = list(contrib = 15), 
                           gradient.cols =  viridis(3,end=0.8,direction = -1),arrowsize=1,labelsize=6,title="")+ 
  theme_light(base_size=24)+
  labs(title = "all lakes")
plot(pcafviz)

pcafviz <- fviz_pca_biplot(pcalakes73cor,axes = c(2,3), repel=T,label = "var",col.ind = "contrib",
                           select.var = list(contrib = 15), 
                           gradient.cols =  viridis(3,end=0.8,direction = -1),arrowsize=1,labelsize=6,title="")+ 
  theme_light(base_size=24)+
  labs(title = "all lakes")

plot(pcafviz)

dev.off()


# -----
# corrplot -----

corlakes73 <- cor(lakes73num, use="complete.obs",method = c("pearson"))
write_xlsx(as.data.frame(corlakes73),"8.version_control/corlakes73.xlsx")
corlakes73_pvalues <- corrplot::cor.mtest(lakes73num,conf.level=0.95)
corrplot::corrplot(corlakes73,method = c("number"),type=c("upper"),tl.cex = 0.6,number.cex=0.45,p.mat = corlakes73_pvalues$p,insig = "blank")

print.cor.signif2(lakes73,"RR")

# -----
# boxplot  -----
lakes73num$trophic <- lakes73$trophic
lakes8num$trophic <- lakes8$trophic
lakes65num$trophic <- lakes65$trophic

# -----
# boxplots lakes 73 ----
pdf("8.version_control/boxplots_lakes73.pdf")

for (i in names(lakes8num)){
  print(ggplot()+geom_boxplot(data = lakes73num, aes(x=i,y=lakes73num[,i]),outlier.shape = NA)+
          geom_jitter(data = lakes65num, aes(x=i,y=lakes65num[,i],col=RR),size=2)+
          geom_jitter(data = lakes8num, aes(x=i,y=lakes8num[,i]),col="red",size=4)+
          ylab(i)+xlab("")+
          scale_color_viridis_c(end = 0.8)+
          theme_light(base_size=20))
}

ggplot(lakes73num)+geom_boxplot(aes(x="RR",y=RR),outlier.shape = NA)+
  geom_jitter(aes(x="RR",y=RR,col=trophic),size=2.5)+
  xlab("")+
  scale_color_viridis_d(end = 0.8)+
  theme_light(base_size=20)

ggplot(lakes73num)+geom_boxplot(aes(x="CNP",y=CNP),outlier.shape = NA)+
  geom_jitter(aes(x="CNP",y=CNP,col=RR),size=2)+
  xlab("")+
  scale_color_viridis_c(end = 0.8)+
  theme_light(base_size=20)

ggplot(lakes73num)+geom_boxplot(aes(x="CP",y=CP),outlier.shape = NA)+
  geom_jitter(aes(x="CP",y=CP,col=RR),size=2)+
  xlab("")+
  scale_color_viridis_c(end = 0.8)+
  theme_light(base_size=20)

ggplot(lakes73num,aes(x="DOC/DP",y=DOC/DP))+geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(col=RR),size=2)+
  xlab("")+
  scale_color_viridis_c(end = 0.8)+
  theme_light(base_size=20)

dev.off()

# -----
# boxplots by trophic state -----
pdf("8.version_control/boxplot_by_trophic.pdf")

for (i in names(lakes8num)){
  print(ggplot()+geom_boxplot(data = lakes73num, aes(x=trophic,y=lakes73num[,i]),outlier.shape = NA)+
          geom_jitter(data = lakes65num, aes(x=trophic,y=lakes65num[,i],col=RR),size=2)+
          geom_jitter(data = lakes8num, aes(x=trophic,y=lakes8num[,i]),col="red",size=4)+
          ylab(i)+xlab("")+
          scale_color_viridis_c(end = 0.8)+
          theme_light(base_size=20))
}

dev.off()

# -----
# other plots -----
ggplot()+geom_text_repel(data=lakes8,aes(x=OD,y=RR*BdgT,col=trophic,label = Lake_ID),nudge_y = 5,size = 5)+
  geom_point(data=lakes8,aes(x=OD,y=RR*BdgT,col=trophic), size = 3)+
  geom_point(data=lakes65,aes(x=OD,y=RR*BdgT,col=trophic),size = 3)+
  theme_light(base_size = 24)

ggplot()+ geom_point(data=lakes65,aes(x=OD,y=RR*BdgT,col=trophic),size = 3)+
  theme_light(base_size = 24)

lm(RR*BdgT~OD,data=lakes73) %>% summary()

ggplot()+geom_point(data=lakes65,aes(x=RR,y=Zn,col="lakes65"))+geom_point(data=lakes8,aes(x=RR,y=Zn,col="outliers"))
ggplot()+geom_point(data=lakes65,aes(x=OD,y=V,col="lakes65"))+geom_point(data=lakes8,aes(x=OD,y=V,col="outliers"))
ggplot()+geom_point(data=lakes65,aes(x=BdgT,y=Cr,col="lakes65"))+geom_point(data=lakes8,aes(x=BdgT,y=Cr,col="outliers"))

ggplot(lakes65)+geom_text(aes(x=lag_bdg,y=RR,label=Lake_ID,col=trophic))


lakes73cor <- lakes73num %>% select(!c("auc","ox_initial","EC_bio","NIVA_date","Long","Altitude","tmax_h","tmax","CNP","Lake_ID","CBA_week","d2H","d18O","Lat","s_275_295","width","DNA","pH_bio","NO2","NO3"))
lakes73cor$RR.DOC <- lakes73cor$RR/lakes73cor$DOC
pdf("8.version_control/comp_cor_RR.pdf",width=10,height = 7)
print.cor.signif2(lakes73cor,"RR",title=F)
print.cor.signif2(lakes73cor,"RR.DOC",title=F)
print.cor.signif2(lakes73cor,"BdgT",title = F)
print.cor.signif2(lakes73cor,"DOC",title=F)
dev.off()




### maps

ggplot(lakes73)+geom_text(aes(x=incub_date,y=RR,col=incub_date,label=Lake_ID),size=5,show.legend = F)+
  theme_light(base_size=20)+theme(legend.position = "none",axis.text.x = element_text(angle = 45))+
  scale_color_brewer(palette="Paired")

ggplot(lakes73)+geom_text(aes(x=incub_date,y=RR,col=DOC,label=Lake_ID),size=5,show.legend = T)+
  theme_light(base_size=20)+theme(axis.text.x = element_text(angle = 30))+
  scale_color_viridis_c(direction = -1, end = 0.8)

ggplot() +
  geom_sf(data=world,fill = "white") +
  coord_sf(xlim = c(4.5,12.9), ylim = c(58, 62), expand = F)+
  geom_point(data=lakes73,aes(x=Long,y=Lat,col=incub_date,size=RR))+
  geom_point(data=lakes8,aes(x=Long,y=Lat),col="black",size=5,shape="*")+
  scale_color_brewer(palette="Paired")+
  theme_void(base_size = 24)+labs(col="Incubation date")

ggplot() +
  geom_sf(data=world,fill = "white") +
  coord_sf(xlim = c(4.5,12.9), ylim = c(58, 62), expand = F)+
  geom_point(data=lakes73,aes(x=Long,y=Lat,col=incub_date,size=DOC))+
  geom_point(data=lakes8,aes(x=Long,y=Lat),col="black",size=5,shape="*")+
  scale_color_brewer(palette="Paired")+
  theme_void(base_size = 24)+labs(col="Incubation date")

ggplot(lakes73num)+geom_text(aes(x=BdgT,y=RR,label=Lake_ID))+theme_light(base_size=20)
ggplot(lakes73num)+geom_text(aes(x=SARuv,y=BdgT,label=Lake_ID))+theme_light(base_size=20)

# -----
# maps -----

pdf("5.100_lakes/RR_maps.pdf",width=10,height=10)

ggplot()+geom_sf(data=norgemap)+
  geom_point(data=lakes73,aes(x=Long,y=Lat,col=RR),size=5)+
  geom_point(data=filter(lakes73,RR >8),aes(x=Long,y=Lat,size=RR))+scale_size(limits=c(8,30),breaks=c(10,20,30),range=c(4,6))+
  theme_void(base_size = 24)+labs(x="",y="",col=expression(atop("RR",paste("(",mu,"mol/L/h)"))),size="outliers")+
  scale_color_gradientn(colors = viridis(9,direction = -1,end=0.9),limits=c(0,8))+
  #geom_text(data=filter(lakes73,RR>8),aes(x=Long,y=Lat,label=round(RR,0)),nudge_y=0.1)+
  guides(color=guide_legend(order=1),size = guide_legend(order=2))+
  theme(legend.position = c(.87,.45))+
  xlim(5,15)+ylim(58,63)

ggplot()+geom_sf(data=norgemap)+
  geom_point(data=lakes73,aes(x=Long,y=Lat,col=RR.OD),size=5)+
  theme_void(base_size = 24)+labs(x="",y="",col=expression(atop("RR/OD",(h^{-1}))))+
  scale_color_viridis_c(direction = -1, end = 0.9)+
  theme(legend.position = c(.9,.45))+
  xlim(5,15)+ylim(58,63)

ggplot()+geom_sf(data=norgemap)+
  geom_point(data=lakes73,aes(x=Long,y=Lat,col=RR.OD),size=5)+
  geom_point(data=filter(lakes73,RR.OD >2),aes(x=Long,y=Lat,size=RR.OD))+scale_size(limits=c(2,10),breaks=c(3,5,10),range=c(4,6))+
  theme_void(base_size = 24)+labs(x="",y="",col=expression(atop("RR/DOC",(h^{-1}))),size="outliers")+
  scale_color_gradientn(colors = viridis(6,direction = -1,end=0.9),limits=c(0,2))+
  #geom_text(data=filter(lakes73,RR.OD >2),aes(x=Long,y=Lat,label=round(RR.OD,2)),nudge_y = -0.09,nudge_x = -0.1)+
  guides(color=guide_legend(order=1),size = guide_legend(order=2))+
  theme(legend.position = c(.87,.45))+
  xlim(5,15)+ylim(58,63)

ggplot()+geom_sf(data=norgemap)+
  geom_point(data=lakes73,aes(x=Long,y=Lat,col=OD.DOC),size=4)+
  theme_void(base_size = 24)+labs(x="",y="",col="%DOC\nconsumed",size="outliers")+
  scale_color_gradientn(colors = viridis(6,direction = -1,end=0.9),limits=c(0,4))+
  geom_point(data=filter(lakes73,OD.DOC >4),aes(x=Long,y=Lat,size=OD.DOC))+scale_size(limits=c(4,20),breaks=c(5,10,15),range=c(4,6))+
  #  geom_text(data=filter(lakes73,OD.DOC >3),aes(x=Long,y=Lat,label=round(OD.DOC,2)),nudge_y = -0.09,nudge_x = -0.1)+  
  guides(color=guide_legend(order=1),size = guide_legend(order=2))+
  theme(legend.position = c(.87,.5))+
  xlim(5,15)+ylim(58,63)

ggplot()+geom_sf(data=norgemap)+
  geom_point(data=lakes73,aes(x=Long,y=Lat,col=OD),size=4)+
  theme_void(base_size = 24)+labs(x="",y="",col="OD\n(umol/L)")+
  scale_color_gradientn(colors = viridis(6,direction = -1,end=0.9),limits=c(0,30))+
  geom_text(data=filter(lakes73,OD >30),aes(x=Long,y=Lat,label=round(OD,0)),nudge_y = -0.09,nudge_x = -0.1)+  
  theme(legend.position = c(.9,.7))+
  xlim(5,15)+ylim(58,63)

ggplot()+geom_sf(data=norgemap)+
  geom_point(data=lakes73,aes(x=Long,y=Lat,col=DN),size=4)+
  geom_point(data=filter(lakes73,DN > 0.6),aes(x=Long,y=Lat,size= DN))+scale_size(limits=c(0.6,0.9))+
  theme_void(base_size = 24)+labs(x="",y="",col="DN\nmg/L",size="Outliers")+
  scale_color_gradientn(colors = viridis(6,direction = -1,end=0.9),limits=c(0,0.6))+
  #  geom_text(data=filter(lakes73,DN > 0.6),aes(x=Long,y=Lat,label=round(DN,1)),nudge_y = -0.09,nudge_x = -0.1)+  
  theme(legend.position = c(.9,.6))+
  xlim(5,15)+ylim(58,63)

ggplot()+geom_sf(data=norgemap)+
  geom_point(data=lakes73,aes(x=Long,y=Lat,col=DP),size=4)+
  geom_point(data=filter(lakes73,DP > 10),aes(x=Long,y=Lat,size= DP))+scale_size(limits=c(10,16))+
  theme_void(base_size = 24)+labs(x="",y="",col="DP\nug/L",size="Outliers")+
  scale_color_gradientn(colors = viridis(6,direction = -1,end=0.9),limits=c(0,10))+
  #  geom_text(data=filter(lakes73,DN > 0.6),aes(x=Long,y=Lat,label=round(DN,1)),nudge_y = -0.09,nudge_x = -0.1)+  
  theme(legend.position = c(.9,.6))+
  xlim(5,15)+ylim(58,63)

ggplot() +
  geom_sf(data=world,fill = "white") +
  coord_sf(xlim = c(4.5,12.9), ylim = c(58, 62), expand = F)+
  geom_point(data=lakes73,aes(x=Long,y=Lat,col=CN),size=3)+
  geom_point(data=filter(lakes73,CN > 0.050),aes(x=Long,y=Lat,size= CN))+scale_size(limits=c(0.05,0.5))+
  scale_color_gradientn(colors = viridis(6,direction = -1,end=0.9),limits=c(0,0.05))+
  theme_void(base_size = 20)+labs(col="CN")


ggplot()+geom_sf(data=norgemap)+
  geom_point(data=lakes73,aes(x=Long,y=Lat,col=CN),size=4)+
  #  geom_point(data=filter(lakes73,CN > 50),aes(x=Long,y=Lat,size= CN))+scale_size(limits=c(50,150),range=c(4,6),breaks=c(50,100,150))+
  theme_void(base_size = 24)+labs(x="",y="",col="CN",size="Outliers")+
  # scale_color_gradientn(colors = viridis(6,direction = -1,end=0.9))+
  scale_color_gradientn(colors = viridis(6,direction = -1,end=0.9),limits=c(8,50))+
  # geom_text(data=filter(lakes73,DN > 0.6),aes(x=Long,y=Lat,label=round(DN,1)),nudge_y = -0.09,nudge_x = -0.1)+ 
  guides(color=guide_legend(order=1),size = guide_legend(order=2))+
  theme(legend.position = c(.9,.6))+
  xlim(5,15)+ylim(58,63)

ggplot()+geom_sf(data=norgemap)+
  geom_point(data=lakes73,aes(x=Long,y=Lat,col=DOC),size=5)+
  #geom_point(data=filter(lakes73,DOC > 20),aes(x=Long,y=Lat,size= DOC))+scale_size(limits=c(20,21))+
  theme_void(base_size = 26)+labs(x="",y="",col="DOC\nmg/L",size="Outliers")+
  # scale_color_gradientn(colors = viridis(6,direction = -1,end=0.9))+
  scale_color_gradientn(colors = viridis(6,direction = -1,end=0.9),limits=c(1,20))+
  #  geom_text(data=filter(lakes73,DOC > 15),aes(x=Long,y=Lat,label=round(DOC,1)),nudge_y = -0.1,nudge_x = -0.0)+ 
  guides(color=guide_legend(order=1),size = guide_legend(order=2))+
  theme(legend.position = c(.75,.2))+
  xlim(5,15)+ylim(58,63)

dev.off()

ggplot() +
  geom_sf(data=world,fill = "white") +
  coord_sf(xlim = c(4.5,12.9), ylim = c(58, 62), expand = F)+
  geom_point(data=lakes73,aes(x=Long,y=Lat,col=Lake_ID),size=3)+
  geom_text_repel(data=lakes73,aes(x=Long,y=Lat,col=Lake_ID,label=Lake_ID),nudge_y=0.1)+
  scale_color_viridis_c()+theme_void(base_size = 20)
ggsave("8.version_control/map_LakeID.jpeg",width=15,height=10)

# -----
# clusters -----

csel <- c("Lake_ID","DOC","TN","TP","a410")
lakes73_clusters <- lakes73[csel] 
lakes73_clusters <- apply(lakes73_clusters,2,as.numeric) %>% as.data.frame() 

clust70c <- lakes73_clusters[,-which(names(lakes73_clusters) %in% c("Lake_ID","RR.DOC","RR","RR.OD","OD","BdgT"))]
rownames(clust70c) <- lakes73_clusters$Lake_ID

dd <- dist(scale(clust70c), method = "euclidian")
hc <- hclust(dd, method = "ward.D2")

dend70c <- as.dendrogram(hc)

# -----
# optimal number of clusters -----
lakes73_scaled <- scale(clust70c) %>% as.data.frame()
lakes73_scaled <- fill.na(lakes73_scaled)

# elbow method
fviz_nbclust(lakes73_scaled, kmeans, method = "wss") 

#silhouette method
fviz_nbclust(lakes73_scaled, kmeans, method = "silhouette")


fviz_nbclust(lakes73_scaled, kmeans, method = "gap_stat")


# -----
# creates clusters ----
nb_clusters <- 4
lakes73$cluster <- cutree(dend70c,nb_clusters)
lakes73_clusters$cluster <- cutree(dend70c,nb_clusters)
colors_cluster <- viridis(nb_clusters,end = 0.8,direction = -1)


dendplotc <- dend70c %>% set("labels_col",value = colors_cluster,k=nb_clusters) %>%
  set("labels_cex",0.8) %>%
  set("branches_k_color", value = colors_cluster,k=nb_clusters) 
plot(dendplotc,horiz = T)

png("8.version_control/dendrogram.png",height = 800, width=300)
plot(dendplotc,horiz = T)
dev.off()


# -----
# graphs and map of clusters -----
pdf("8.version_control/bdg_by_cluster.pdf",width=10,height=8)
ggplot(lakes73,aes(y=RR.OD,x=trophic))+geom_boxplot(aes(group=trophic),outlier.shape = NA)+geom_jitter(aes(col=DOC),width=0.3,size=3)+
  scale_color_viridis_c(direction = -1,end=0.9)+theme_light(base_size=26)+ylab("RR/OD")
ggplot(lakes73,aes(y=RR,x=trophic))+geom_boxplot(aes(group=trophic),outlier.shape = NA)+geom_jitter(aes(col=DOC),width=0.3,size=3)+
  scale_color_viridis_c(direction = -1,end=0.9)+theme_light(base_size=26)+ylab("RR")
ggplot(lakes73,aes(y=RR.DOC,x=trophic))+geom_boxplot(aes(group=trophic),outlier.shape = NA)+geom_jitter(aes(col=DOC),width=0.3,size=3)+
  scale_color_viridis_c(direction = -1,end=0.9)+theme_light(base_size=26)+ylab("RR/DOC")

dev.off()

ggplot()+geom_sf(data=norgemap)+
  geom_point(data=lakes73,aes(x=Long,y=Lat,col=as.character(cluster)),size=6)+
  theme_void(base_size = 26)+labs(x="",y="",col="cluster")+
  scale_color_manual(values=colors_cluster)+
  theme(legend.position = c(.75,.15))+
  xlim(5,15)+ylim(58,63)
ggsave("8.version_control/map_clusters.png",width=10,height = 10)


# -----
# ANOVA -----

a1 <- aov(RR.OD~as.factor(cluster),data=lakes73)
posthoc <- TukeyHSD(x=a1,conf.level=0.95)

aov_summary <- as.data.frame(names(lakes73_clusters))
aov_summary$pvalue <- NA
aov_summary$effect_size <- NA
aov_summary$ptt12 <- NA
aov_summary$ptt13 <- NA
aov_summary$ptt23 <- NA

for (i in 1:length(names(lakes73_clusters))){
  param <- names(lakes73_clusters)[i]
  aovt <- aov(lakes73_clusters[,param]~lakes73_clusters$cluster)
  pval <- summary(aovt)[[1]]$`Pr(>F)`[1]
  es <- eta_squared(aovt)[2]
  ptt <- pairwise.t.test(lakes73_clusters[,param],lakes73_clusters$cluster,p.adjust.method = "bonferroni")
  aov_summary$pvalue[i] <- pval
  aov_summary$effect_size[i] <- es[1,1]
  aov_summary$ptt12[i] <- ptt$p.value[1]
  aov_summary$ptt13[i] <- ptt$p.value[2]
  aov_summary$ptt23[i] <- ptt$p.value[4]
}

names(aov_summary)[1] <- "param"
write_xlsx(aov_summary,path="8.version_control/anovas_clusters_fdr.xlsx")


ptt <- aov_summary %>% select(c("param","ptt12","ptt13","ptt23","pvalue"))

ptt$pvalue[which(ptt$pvalue > 0.05)] <- NA
ptt$pvalue[which(ptt$pvalue < 0.05)] <- "significant"
ptt$pvalue[which(is.na(ptt$pvalue))] <- "non-significant"

ptt$ptt12[which(ptt$ptt12 > 0.05)] <- NA
ptt$ptt12[which(ptt$ptt12 < 0.05)] <- "Cluster 1 vs Cluster 2"

ptt$ptt13[which(ptt$ptt13 > 0.05)] <- NA
ptt$ptt13[which(ptt$ptt13 < 0.05)] <- "Cluster 1 vs Cluster 3"

ptt$ptt23[which(ptt$ptt23 > 0.05)] <- NA
ptt$ptt23[which(ptt$ptt23 < 0.05)] <- "Cluster 2 vs Cluster 3"


ptt.graph <- ptt %>% filter(!param %in% c("Sample_ID","cluster"))


x <- c("RR.OD", "CN","pH","Al",
       "TOC","DOC","Fe","Altitude","Co59","V51","Temperature","Cu63",
       "Na","DN","Mg","K","Ni58",
       "Alkalinity","Conductivity","CO2","TN","TP","Cl-","SO4","BdgT",
       "SUVA", "As75",
       "CH4","cells_counts.mL","SR","CP",
       "DP","Ca","O2","tmax","RR")
x.graph <- x[x %in% ptt.graph$param]

ggplot(ptt.graph,aes(y=param))+geom_point(aes(x=ptt12,col="12",shape=pvalue),size=5)+
  geom_point(aes(x=ptt13,col="13",shape=pvalue),size=5)+
  geom_point(aes(x=ptt23,col="23",shape=pvalue),size=5)+scale_x_discrete(position = "top")+
  labs(x="",y="",col="",shape = "ANOVA p-value")+scale_color_viridis_d(end=0.8)+theme_light(base_size = 24)+
  theme(panel.grid.major.x = element_blank(),legend.position = "bottom",axis.text.x = element_text(angle = 20,hjust=0,vjust=1),
        plot.margin = unit(c(0,7,0,0),"cm"))+
  scale_y_discrete(limits = (rev(x.graph)))
ggsave("5. 100 lakes/cluster_difference.png",height = 15,width = 10)

# -----
# Kruskal Wallis -----

k1 <- kruskal.test(DOC~cluster,data=lakes73)
pwt <- pairwise.wilcox.test(lakes73$RR,lakes73$cluster,p.adjust.method = "holm")

kw_summary <- as.data.frame(names(lakes73))
kw_summary$pvalue <- NA
#aov_summary$effect_size <- NA
kw_summary$pwt12 <- NA
kw_summary$pwt13 <- NA
kw_summary$pwt23 <- NA

for (i in 1:length(names(lakes73))){
  param <- names(lakes73)[i]
  kwt <- kruskal.test(lakes73[,param]~lakes73$cluster)
  pval <- kwt$p.value
  #es <- eta_squared(aovt)$Eta_Sq_partial
  pwt <- pairwise.wilcox.test(lakes73[,param],lakes73$cluster,p.adjust.method = "holm",paired=FALSE)
  kw_summary$pvalue[i] <- pval
  #kw_summary$effect_size[i] <- es
  kw_summary$pwt12[i] <- pwt$p.value[1]
  kw_summary$pwt13[i] <- pwt$p.value[2]
  kw_summary$pwt23[i] <- pwt$p.value[4]
}

names(kw_summary)[1] <- "param"
write_xlsx(kw_summary,path="5. 100 lakes/kw_clusters_fdr.xlsx")


pwt <- kw_summary %>% select(c("param","pwt12","pwt13","pwt23","pvalue"))

pwt$pvalue[which(pwt$pvalue > 0.05)] <- NA
pwt$pvalue[which(pwt$pvalue < 0.05)] <- "significant"
pwt$pvalue[which(is.na(pwt$pvalue))] <- "non-significant"

pwt$pwt12[which(pwt$pwt12 > 0.05)] <- NA
pwt$pwt12[which(pwt$pwt12 < 0.05)] <- "Cluster 1 ≠ Cluster 2"

pwt$pwt13[which(pwt$pwt13 > 0.05)] <- NA
pwt$pwt13[which(pwt$pwt13 < 0.05)] <- "Cluster 1 ≠ Cluster 3"

pwt$pwt23[which(pwt$pwt23 > 0.05)] <- NA
pwt$pwt23[which(pwt$pwt23 < 0.05)] <- "Cluster 2 ≠ Cluster 3"


pwt.graph <- pwt %>% filter(!param %in% c("Sample_ID","cluster"))

# -----
# graph -----
x <- c("temperature","Al","SUVA","Fe","CH4","Cl","Co","Cu",
       "RR.OD", "CN","CP","pH","SR","BdgT",
       "TOC","DOC","RR.OD","OD.DOC","Altitude",
       "Na","V","As","Ni",
       "DN","Mg","K",
       "Alkalinity","Conductivity","CO2","TN","TP","SO4","DP",
      "Ca","cells",
      "O2","tmax","RR",
      "dist_bact")
x.graph <- x[x %in% pwt.graph$param]

ggplot(pwt.graph,aes(y=param))+geom_point(aes(x=pwt12,col="12",shape=pvalue),size=5)+
  geom_point(aes(x=pwt13,col="13",shape=pvalue),size=5)+
  geom_point(aes(x=pwt23,col="23",shape=pvalue),size=5)+scale_x_discrete(position = "top")+
  labs(x="",y="",col="",shape = "Kruskal-Wallis p-value")+scale_color_viridis_d(end=0.8)+theme_light(base_size = 24)+
  theme(panel.grid.major.x = element_blank(),legend.position = "bottom",axis.text.x = element_text(angle = 20,hjust=0,vjust=1),
        plot.margin = unit(c(0,7,0,0),"cm"))+
  scale_y_discrete(limits = (rev(x.graph)))
ggsave("5. 100 lakes/cluster_difference.png",height = 15,width = 10)

# -----
# data presentation ----
library(grid)
vp <- viewport(x=0.45,y=0.5,width=0.9,height=0.9)

pdf("5. 100 lakes/CBA_figures.pdf",width=10,height = 5)
pushViewport(vp)
A <- ggplot(lakes73)+geom_boxplot(aes(x=RR,y="RR\n(uM/h)"))+theme_light(base_size=15)+labs(x="",y="")
A2 <- ggplot(lakes73)+geom_boxplot(aes(x=OD,y="OD\n(uM)"))+theme_light(base_size=15)+labs(x="",y="")
B <- ggplot(lakes73)+geom_boxplot(aes(x=RR.OD,y="RR_s\n(h-1)"))+theme_light(base_size=15)+labs(x="",y="")
C <- ggplot(lakes73)+geom_boxplot(aes(x=RR.OD,y="RR_OD\n(h-1)"))+theme_light(base_size=15)+labs(x="",y="")
gA <- ggplotGrob(A)
gA2 <- ggplotGrob(A2)
gB <- ggplotGrob(B)
gC <- ggplotGrob(C)
grid::grid.newpage()
grid::grid.draw(cbind(rbind(gA,gA2),rbind(gB,gC)))

D <- ggplot(lakes73)+geom_boxplot(aes(x=DOC,y="DOC (mg/L)"))+theme_light(base_size=15)+labs(x="",y="")
E <- ggplot(lakes73)+geom_boxplot(aes(x=DN,y="DN (mg/L)"))+theme_light(base_size=15)+labs(x="",y="")
F <- ggplot(lakes73)+geom_boxplot(aes(x=DP,y="DP (ug/L)"))+theme_light(base_size=15)+labs(x="",y="")
G <- ggplot(lakes73)+geom_boxplot(aes(x=TOC,y="TOC (mg/L)"))+theme_light(base_size=15)+labs(x="",y="")
H <- ggplot(lakes73)+geom_boxplot(aes(x=TN,y="TN (mg/L)"))+theme_light(base_size=15)+labs(x="",y="")
I <- ggplot(lakes73)+geom_boxplot(aes(x=TP,y="TP (ug/L)"))+theme_light(base_size=15)+labs(x="",y="")
J <- ggplot(lakes73)+geom_boxplot(aes(x=DOC/DN,y="C:N (dissolved)"))+theme_light(base_size=15)+labs(x="",y="")
K <- ggplot(lakes73)+geom_boxplot(aes(x=TOC/TN,y="C:N (total)"))+theme_light(base_size=15)+labs(x="",y="")
gD <- ggplotGrob(D)
gE <- ggplotGrob(E)
gF <- ggplotGrob(F)
gG <- ggplotGrob(G)
gH <- ggplotGrob(H)
gI <- ggplotGrob(I)
gJ <- ggplotGrob(J)
gK <- ggplotGrob(K)
grid::grid.newpage()
grid::grid.draw(cbind(rbind(gD,gE,gF,gJ),rbind(gG,gH,gI,gK)))

L <- ggplot(lakes73)+geom_boxplot(aes(x=pH_lab,y="pH"))+theme_light(base_size=15)+labs(x="",y="") 
gL <- ggplotGrob(L)
M <- ggplot(lakes73)+geom_boxplot(aes(x=alk_meq.L,y="alkalinity (meq/L)"))+theme_light(base_size=15)+labs(x="",y="")
gM <- ggplotGrob(M)
N <- ggplot(lakes73)+geom_boxplot(aes(x=Altitude,y="Altitude (m)"))+theme_light(base_size=15)+labs(x="",y="")
grid::grid.newpage()
grid::grid.draw(rbind(gL,gM))

N <- ggplot(lakes73)+geom_boxplot(aes(x=Altitude,y="Altitude (m)"))+theme_light(base_size=15)+labs(x="",y="")
gN <- ggplotGrob(N)

grid::grid.newpage()
grid::grid.draw(cbind(rbind(gD,gE,gJ),rbind(gF,gL,gN)))


ggplot(lakes73)+geom_point(aes(x=OD,y=RR,col=DOC))+scale_color_viridis_c(direction = -1, end = 0.9)+theme_light(base_size=20)

dev.off()

summary_lakes73 <- do.call(cbind,lapply(lakes73,summary)) %>% as.data.frame()
summary_lakes73[,1] <- rownames(summary_lakes73)
write_xlsx(summary_lakes73,"5. 100 lakes/summary_lakes73.xlsx")


# -----
# correLations ----
pdf("5.100_lakes/all.parameters.correLations.pdf",width=15,height = 8)
sapply(names(lakes73),print.cor.signif2, df=lakes73)
dev.off()

# -----
# print cor all lakes -----

pdf("8.version_control/corplots.pdf",width=10,height = 8)

lakes73num$trophic <- NULL

#print.cor.signif2(lakes73num,"RR")
print.cor.signif2(select(lakes73num,c("RR","TOC","CN","CP","DOC","tmax","EC_bio","ox_initial","OD")),"RR",title =  "73 lakes")

#print.cor.signif2(lakes73num,"OD")
print.cor.signif2(select(lakes73num,c("RR","TOC","CN","CP","DOC","tmax","OD")),"OD",title = "73 lakes")

#print.cor.signif(lakes73num,"BdgT")
print.cor.signif2(select(lakes73num,c("BdgT","DN","TN","NO2","NO3","K","Ca","pH","Alkalinity","EC_Kje","c_O2","Mg","TP","Cd")),"BdgT",title = "73 lakes")

# -----
# print cor without outliers -----

lakes65num <- lakes65 %>% sapply(as.numeric) %>% as.data.frame() %>% select_if(not_all_na)

#print.cor.signif2(lakes65num,"RR")
print.cor.signif2(select(lakes65num,c("RR","OD","CN","K","Na","EC","Mn","Mg","As")),"RR",title = "Without outliers")
#print.cor.signif(lakes65num,"OD",title = F)
print.cor.signif2(select(lakes65num,c("OD","RR","BdgT","c_O2","pH","tmax","Co")),"OD",title = "Without outliers")
#print.cor.signif(lakes65num,"BdgT")
print.cor.signif2(select(lakes65num,c("BdgT","OD","DN","NO2","TN","K","NO3","Ca","c_O2","EC","pH","Mg","Alkalinity","TP","Hg_total","As","CN","Cells")),"BdgT",title = "Without outliers")

# -----
# print cor outliers -----
lakes8 <- lakes73 %>% filter(RR > 10)
lakes8num <- lakes8 %>% sapply(as.numeric) %>% as.data.frame() %>% select_if(not_any_na)

#print.cor.signif(lakes8num,"RR")
print.cor.signif2(select(lakes8num,c("RR","OD","SARuv","Fe","Zn","Pb")),"RR",title = "Outliers")
#print.cor.signif(lakes8num,"OD")
print.cor.signif2(select(lakes8num,c("OD","tmax","TOC","RR","CP","BdgT","V")),"OD",title = "Outliers")
#print.cor.signif(lakes8num,"BdgT")
print.cor.signif2(select(lakes8num,c("BdgT","Cr","TOC","CP","DOC","CN","OD")),"BdgT",title = "Outliers")

dev.off()
# -----
# print cor all lakes spearman -----

pdf("8.version_control/corplots_spearman.pdf",width=10,height = 8)

lakes73num$trophic <- NULL

#print.cor.signif2(lakes73num,"RR",meth="spearman")
print.cor.signif2(select(lakes73num,c("RR","OD","CN","EC","Cl","K","SO4","Hg_total","Mn","NO3","Na","Mg")),"RR",title =  "73 lakes", meth = "spearman")

#print.cor.signif2(lakes73num,"OD",meth = "spearman")
print.cor.signif2(select(lakes73num,c("RR","BdgT","c_O2","Na","Hg_total","p_N2","OD")),"OD",title = "73 lakes", meth = "spearman")

#print.cor.signif2(lakes73num,"BdgT", meth = "spearman")
print.cor.signif2(select(lakes73num,c("BdgT","pH","Alkalinity","OD","Ca","O2","SO4","NO2","Fe","CN","CP","SARvis")),"BdgT",title = "73 lakes", meth = "spearman")

# -----
# print cor without outliers spearman  -----

lakes65num <- lakes65 %>% sapply(as.numeric) %>% as.data.frame() %>% select_if(not_all_na)

#print.cor.signif2(lakes65num,"RR",meth="spearman")
print.cor.signif2(select(lakes65num,c("RR","OD","CN","K","Na","EC","Mn","Mg","NO3","Hg_total")),"RR",title = "Without outliers", meth = "spearman")
#print.cor.signif2(lakes65num,"OD",title = "without ouliers",meth = "spearman")
print.cor.signif2(select(lakes65num,c("OD","RR","BdgT","pH","O2","SARuv","Hg_total")),"OD",title = "Without outliers", meth = "spearman")
#print.cor.signif2(lakes65num,"BdgT", meth = "spearman")
print.cor.signif2(select(lakes65num,c("BdgT","pH","OD","Alkalinity","Ca","O2","SO4","NO2","Mg","DP","Cells","Fe","CP","CN","SARvis","EC_bio")),"BdgT",title = "Without outliers", meth = "spearman")

# -----
# print cor outliers spearman -----
lakes8 <- lakes73 %>% filter(RR > 10)
lakes8num <- lakes8 %>% sapply(as.numeric) %>% as.data.frame() %>% select_if(not_any_na)

#print.cor.signif2(lakes8num,"RR", meth = "spearman")
print.cor.signif2(select(lakes8num,c("RR","SARuv","V","Zn","Co")),"RR",title = "Outliers", meth = "spearman")
#print.cor.signif2(lakes8num,"OD", meth = "spearman")
print.cor.signif2(select(lakes8num,c("OD","CP","TOC","BdgT","CP","BdgT","TP")),"OD",title = "Outliers", meth = "spearman")
#print.cor.signif2(lakes8num,"BdgT",meth = "spearman")
print.cor.signif2(select(lakes8num,c("BdgT","p_CO2","OD","c_CO2")),"BdgT",title = "Outliers", meth = "spearman")

dev.off()


# -----
# print correlation by trophic state -----
lakes73num$trophic <- lakes73$trophic

pdf("8.version_control/cor_trophic_state.pdf",width=10,height=8)

lakes73eu <- lakes73num %>% filter(trophic == "eu") %>% select_if(not_any_na)
lakes73eu$trophic <- NULL
sapply(names(lakes73eu),print.cor.signif,df = lakes73eu, title = F)


lakes73oligo <- lakes73num %>% filter(trophic == "oligo") %>% select_if(not_any_na)
lakes73oligo$trophic <- NULL
sapply(names(lakes73oligo),print.cor.signif,df = lakes73oligo, title = F)


lakes73dys <- lakes73num %>% filter(trophic == "dys") %>% select_if(not_any_na)
lakes73dys$trophic <- NULL
sapply(names(lakes73dys),print.cor.signif,df = lakes73dys, title = F)

dev.off()

pdf("8.version_control/corbytrophicstate.pdf",width=10,height=8)
print.cor.signif2(select(lakes73eu,c("RR","OD","CN","CP","DOC","SR","TOC","s_275_295","SUVA","DP","V")),"RR",title = F)
print.cor.signif2(select(lakes73eu,c("RR","OD","CN","CP","DOC","SR","TOC","Cr","SR","SUVA","V")),"OD",title = F)
print.cor.signif2(select(lakes73eu,c("BdgT","Zn","TP","DP","DN","pH_bio")),"BdgT",title = F)

print.cor.signif2(select(lakes73oligo,c("OD","RR","CN","H","DOC","CP","pH")),"RR",title = F)
print.cor.signif2(select(lakes73oligo,c("OD","RR","CN","H","DOC","CP","c_O2")),"OD",title = F)
print.cor.signif2(select(lakes73oligo,c("BdgT","Ca","alkalinity","Mn","DN","SARuv","NO2","conductivity","TN","SO4","c_O2","pH")),"BdgT",title = F)

print.cor.signif2(select(lakes73dys,c("RR","OD","TOC","CN","CP")),"RR",title = F)
print.cor.signif2(select(lakes73dys,c("RR","OD","TOC","CN","CP")),"OD",title = F)
print.cor.signif2(select(lakes73dys,c("BdgT","Ca","conductivity","Na","SO4","NO2","B","K","Cl","Mg","NO3","pH","TN","DN","c_CO2","Cu","As","CN","NP")),"BdgT",title = F)
dev.off()

# -----
# correlations - figures ----

pdf("5. 100 lakes/RR_corr.pdf",height = 10,width=15)
cor_param_RR <- c("Sample_ID","OD","RR","CN","CP") 
cor_RR <- lakes73[cor_param_RR] %>% setNames(c("Sample_ID","OD","RR","C:N","C:P"))
print.cor.signif2(cor_RR,"RR",title=F,y1=-0.5)

cor_param_RR.OD <- c("Sample_ID","OD","RR","s_275_295","SR","DN","Al_NIVA","SUVA","DOC","RR.OD","OD.DOC")
cor_RR.OD <- lakes73[cor_param_RR.OD] %>% setNames(c("Sample_ID","OD","RR","s_275_295","SR","DN","Al","SUVA","DOC","RR/DOC","%DOM"))
print.cor.signif2(cor_RR.OD,"RR/DOC",title=F,y1=-0.5)
print.cor.signif2(cor_RR.OD,"%DOM",title=F,y1=-0.5,y2=1.3)

cor_param_RR.OD <- c("Sample_ID","Al_NIVA","Mg","K","DN","alkalinity","pH","RR.OD","CN")
cor_RR.OD <- lakes73[cor_param_RR.OD] %>% setNames(c("Sample_ID","Al","Mg","K","DN","alkalinity","pH","RR/OD","C:N"))
print.cor.signif2(cor_RR.OD,"RR/OD",title = F,y1=-0.7,y2=0.5)

cor_param_DOC <- c("DOC","TOC","V","Al","SUVA","Cu","Na","Fe","Co","Pb","Zn","cells","Br","DN","TN","CP","CN","K","conductivity","Ni","Mg","CH4","RR.OD","OD.DOC","temperature","SR","pH","Altitude")
cor_DOC <- lakes73[cor_param_DOC]
print.cor.signif2(cor_DOC,"DOC",title = F)

dev.off()


# -----
# linear regression ----

lm(RR~TN*TP*TOC*Ca*pH,data=lakes73) %>% summary()
lm(RR~CN*CP,data=lakes73) %>% summary()
lm(RR~CN*CP+pH,data=lakes73) %>% summary()
lm(RR~CN*CP*pH,data=lakes73) %>% summary()
lm(RR~CN*CP*Ca,data=lakes73) %>% summary()
lm(RR~CN*CP*Ca*pH*a600,data=lakes73) %>% summary()
lm(RR~CN*CP*Ca*pH,data=lakes73) %>% summary()
lm(RR~CN*CP*pH,data=lakes65) %>% summary()
lm(RR~CN*CP*pH,data=lakes8) %>% summary()



# -----
# Figures for paper -----
# -----
# PCA -----
pca_params <- c("pH","SR","RR.OD","DN","SUVA","DOC","alkalinity","RR","RR.OD","CN","CP")
lakes73pca <- lakes73[pca_params]
pca70lakes <- prcomp(~.,data=lakes73pca,scale.=T,center=T)
pcafviz <- fviz_pca_biplot(pca70lakes, repel=T,label = "var",col.var = "contrib",col.ind = "contrib",
                           gradient.cols =  viridis(3,end=0.8,direction = -1),arrowsize=1,labelsize=6,title="")+ theme_light(base_size=24) 
pcafviz23 <- fviz_pca_biplot(pca70lakes, axes = c(2,3),repel=T,label = "var",col.var = "contrib",col.ind = "contrib",
                             gradient.cols =  viridis(3,end=0.8,direction = -1),arrowsize=1,labelsize=6,title="")+ theme_light(base_size=24) 
plot(pcafviz)

ggsave("5. 100 lakes/RR_pca.png",height = 8, width = 10)
plot(pcafviz23)

pca70lakes.rot <- as.data.frame(summary(pca70lakes)$rotation)
pca70lakes.rot$param <- rownames(pca70lakes.rot)

PC1df <- pca70lakes.rot %>% select(c("param","PC1")) %>% arrange(desc(PC1))
PC2df <- pca70lakes.rot %>% select(c("param","PC2")) %>% arrange(desc(PC2))
PC3df <- pca70lakes.rot %>% select(c("param","PC3")) %>% arrange(desc(PC3))

gPC1 <- ggplot(PC1df,aes(y=reorder(param,order(PC1,decreasing=F)),x=abs(PC1),label=round(PC1,2)))+
  geom_col(aes(fill=PC1),show.legend = F)+
  geom_text(nudge_x = 0.04)+
  theme_bw(base_size=20)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_fill_gradient(high="red",low="skyblue")+labs(x="",y="PC1 (31%)")
gPC2 <- ggplot(PC2df,aes(y=reorder(param,order(PC2,decreasing=F)),x=abs(PC2),label=round(PC2,2)))+
  geom_col(aes(fill=PC2),show.legend = F)+ 
  geom_text(nudge_x = 0.04)+
  theme_bw(base_size=20)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_fill_gradient(high="red",low="skyblue")+labs(x="",y="PC2 (21%)")#+
gPC3 <- ggplot(PC3df,aes(y=reorder(param,order(PC3,decreasing=F)),x=abs(PC3),label=round(PC3,2)))+
  geom_col(aes(fill=PC3),show.legend = F)+ 
  geom_text(nudge_x = 0.04)+
  theme_bw(base_size=20)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_fill_gradient(high="red",low="skyblue")+labs(x="",y="PC3 (11%)")


pdf("5. 100 lakes/pca_postcluster.pdf",width=12,heigh=8)
plot(pcafviz)
plot(pcafviz23)
plot(gPC1)
plot(gPC2)
plot(gPC3)
dev.off()












# -----
# extra plots ----


ggplot(lakes73,aes(x=DOC,y=RR.OD))+geom_point()+geom_smooth(method = "gam",formula = y~s(x))+theme_light(base_size=24)
ggplot(lakes73,aes(x=DN,y=RR.OD))+geom_point()+geom_smooth(method = "gam",formula = y~s(x))+theme_light(base_size=24)
ggplot(lakes73,aes(x=CN,y=RR.OD))+geom_point()+geom_smooth(method = "gam",formula = y~s(x))+theme_light(base_size=24)
ggplot(lakes73,aes(x=pH,y=RR.OD))+geom_point()+geom_smooth(method = "gam",formula = y~s(x))+theme_light(base_size=24)


ggplot(lakes73,aes(x=DOC,y=RR))+geom_point()+geom_smooth(method = "gam",formula = y~s(x))+theme_light(base_size=24)
ggplot(lakes73,aes(x=DN,y=RR))+geom_point()+geom_smooth(method = "gam",formula = y~s(x))+theme_light(base_size=24)
ggplot(lakes73,aes(x=CN,y=RR))+geom_point()+geom_smooth(method = "gam",formula = y~s(x))+theme_light(base_size=24)


# ---
# lasso model -----
# install.packages("glmnet", repos = "http://cran.us.r-project.org")
library("glmnet")
library(glmnetUtils)
library(mice)
library(broom.mixed)

# create multiple imputation datasets to fill NA -----
param2remove <- c("auc","width","OD","Lake_ID","NVE_number","NIVA_date","X","Long","Lat","Altitude","CBA_week","tmax","tmax_h","ox_initial","F","pH_bio","EC_bio","NO2","NO3","DNA","p_N2","p_O2","c_O2","p_CO2","c_CO2","p_CH4","p_N2O","d18O","d2H","Tso","a254","a400","a410","a600","a275","a295","s_275_295","a350","a400","s_350_400","lag_bdg","H","NP","CNP","dist_bact")
lakes73clean <- lakes73num %>% dplyr::select(!which(names(lakes73num)%in% param2remove))
write_xlsx(lakes73clean,"8.version_control/lakes73clean.xlsx")
M <- 10
RRmice <- lakes73clean %>% mice(method = "cart",m = M)
# -----
# cor.test on mice -----
predvar <- names(lakes73clean)[-which(names(lakes73clean)=="RR")]

pearson_pooled <- predvar %>% as.data.frame() %>% setNames("param")
estimates <- c()
p_values <- c()
for(y in 1:length(predvar)){
  cor_formula <- as.formula(paste("~RR+",predvar[y],sep=""))
  for(x in 1:M){
  df <- complete(RRmice,x)
  dfcor <- cor.test(formula = cor_formula, data = df)
  estimates[x] <- dfcor$estimate[[1]]
  p_values[x] <- dfcor$p.value[[1]]
}

cor_pooled <- pool.scalar(estimates,p_values,n=73,k=1) 
pearson_pooled$r[y] <- cor_pooled$qbar
pearson_pooled$p_value[y] <- cor_pooled$ubar

}

g <- ggplot(data = filter(pearson_pooled,p_value <= 0.05),aes(x=reorder(param,order(r,decreasing=F)),y=r,label=param))+
  geom_col(aes(fill=p_value))+ 
  geom_text(angle=90,hjust=0,nudge_y = 0.01,size=7)+
  theme_bw(base_size=28)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.grid.minor = element_blank())+
  scale_fill_gradientn(colors = c("skyblue","dodgerblue4","firebrick"),
                       breaks = c(0,0.01,0.025, 0.05),labels = c("0","0.01","0.025","0.05"),limits = c(0,0.05), 
                       guide="legend")+
  labs(x = "", y= "r2", title = "Pearson correlation coefficient with RR") + ylim(y1,y2)
  print(g)


# -----
# cor.test on mice for BdgT -----

predvar <- names(lakes73clean)[-which(names(lakes73clean)=="BdgT")]

pearson_pooled_bdgT <- predvar %>% as.data.frame() %>% setNames("param")
estimates <- c()
p_values <- c()
for(y in 1:length(predvar)){
  cor_formula <- as.formula(paste("~BdgT+",predvar[y],sep=""))
  for(x in 1:M){
    df <- filter(RRmice_long,.imp == x)
    dfcor <- cor.test(formula = cor_formula, data = df)
    estimates[x] <- dfcor$estimate[[1]]
    p_values[x] <- dfcor$p.value[[1]]
  }
  
  cor_pooled <- pool.scalar(estimates,p_values,n=73,k=1) 
  pearson_pooled_bdgT$r[y] <- cor_pooled$qbar
  pearson_pooled_bdgT$p_value[y] <- cor_pooled$ubar
  
}

g <- ggplot(data = filter(pearson_pooled_bdgT,p_value < 0.05),aes(x=reorder(param,order(r,decreasing=T)),y=r,label=param))+
  geom_col(aes(fill=p_value))+ 
  geom_text(angle=90,hjust=0,nudge_y = 0.01,size=7)+
  theme_bw(base_size=28)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.grid.minor = element_blank())+
  scale_fill_gradientn(colors = c("skyblue","dodgerblue4","firebrick"),
                       breaks = c(0,0.01,0.025, 0.05),labels = c("0","0.01","0.025","0.05"),limits = c(0,0.05), 
                       guide="legend")+
  xlab("")+ ylim(y1,y2)
print(g)
# -----
# cor.test on mice with Spearman -----
predvar <- names(lakes73clean)[-which(names(lakes73clean)=="RR")]

pearson_pooled <- predvar %>% as.data.frame() %>% setNames("param")
estimates <- c()
p_values <- c()
for(y in 1:length(predvar)){
  cor_formula <- as.formula(paste("~RR+",predvar[y],sep=""))
  for(x in 1:M){
    df <- complete(RRmice,x)
    dfcor <- cor.test(formula = cor_formula, data = df, method = "spearman")
    estimates[x] <- dfcor$estimate[[1]]
    p_values[x] <- dfcor$p.value[[1]]
  }
  
  cor_pooled <- pool.scalar(estimates,p_values,n=73,k=1) 
  pearson_pooled$r[y] <- cor_pooled$qbar
  pearson_pooled$p_value[y] <- cor_pooled$ubar
  
}

g <- ggplot(data = filter(pearson_pooled,p_value < 0.05),aes(x=reorder(param,order(r,decreasing=T)),y=r,label=param))+
  geom_col(aes(fill=p_value))+ 
  geom_text(angle=90,hjust=0,nudge_y = 0.01,size=7)+
  theme_bw(base_size=28)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.grid.minor = element_blank())+
  scale_fill_gradientn(colors = c("skyblue","dodgerblue4","firebrick"),
                       breaks = c(0,0.01,0.025, 0.05),labels = c("0","0.01","0.025","0.05"),limits = c(0,0.05), 
                       guide="legend")+
  xlab("")+ ylim(y1,y2)
print(g)

df2 <- filter(pearson_pooled,p_value < 0.05)
df2$param <- as.factor(df2$param)
g <- ggplot(data = df2,aes(x=reorder(param,order(df2$r,decreasing = T)),y=r,label=param))+
  geom_col(aes(fill=p_value))+ 
  geom_text(angle=90,hjust=0,nudge_y = 0.01,size=7)+
  theme_bw(base_size=28)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.grid.minor = element_blank())+
  scale_fill_gradientn(colors = c("skyblue","dodgerblue4","firebrick"),
                       breaks = c(0,0.01,0.025, 0.05),labels = c("0","0.01","0.025","0.05"),limits = c(0,0.05), 
                       guide="legend")+
  xlab("")+ ylim(y1,y2)
print(g)


# -----
# cor.test on mice for BdgT -----

predvar <- names(lakes73clean)[-which(names(lakes73clean)=="BdgT")]

pearson_pooled_bdgT <- predvar %>% as.data.frame() %>% setNames("param")
estimates <- c()
p_values <- c()
for(y in 1:length(predvar)){
  cor_formula <- as.formula(paste("~BdgT+",predvar[y],sep=""))
  for(x in 1:M){
    df <- filter(RRmice_long,.imp == x)
    dfcor <- cor.test(formula = cor_formula, data = df)
    estimates[x] <- dfcor$estimate[[1]]
    p_values[x] <- dfcor$p.value[[1]]
  }
  
  cor_pooled <- pool.scalar(estimates,p_values,n=73,k=1) 
  pearson_pooled_bdgT$r[y] <- cor_pooled$qbar
  pearson_pooled_bdgT$p_value[y] <- cor_pooled$ubar
  
}

g <- ggplot(data = filter(pearson_pooled_bdgT,p_value < 0.05),aes(x=reorder(param,order(r,decreasing=T)),y=r,label=param))+
  geom_col(aes(fill=p_value))+ 
  geom_text(angle=90,hjust=0,nudge_y = 0.01,size=7)+
  theme_bw(base_size=28)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.grid.minor = element_blank())+
  scale_fill_gradientn(colors = c("skyblue","dodgerblue4","firebrick"),
                       breaks = c(0,0.01,0.025, 0.05),labels = c("0","0.01","0.025","0.05"),limits = c(0,0.05), 
                       guide="legend")+
  xlab("")+ ylim(y1,y2)
print(g)

# -----
# mice and lasso -----
# test 
set.seed(1)
test <- cv.glmnet(formula = RR~.,data=lakes73clean)
test2 <- glmnet(formula = RR~.,data=lakes73clean)


# lasso for RR -----
lasso_coef <- c("Intercept",names(dplyr::select(lakes73clean,!c("RR")))) %>% as.data.frame() %>% setNames("param")

for(z in 1:M){
  lasso_coef %>% tibble::add_column(z = NA)
  data <- complete(RRmice,z)
  set.seed(1)
  fit <- cv.glmnet(formula = RR~., data = data)
  coef_fit <- fit %>% coef(s="lambda.min") %>% summary()
  dev_fit <- fit$glmnet.fit %>% deviance()
  for (w in 1:dim(coef_fit)[1]){
    lasso_coef[coef_fit[w,1],z+1] <- coef_fit[w,3]
  }
}

lasso_coef$pooled <- NA
for (x in 1:length(lasso_coef$param)){
  estimates <- lasso_coef[x,] %>% dplyr::select(c(2:M+1)) %>% as.numeric()
  lasso_coef$pooled[x] <- pool.scalar(estimates,c(rep(1,10)),n=73,k=1)$qbar
}

write_xlsx(lasso_coef,"8.version_control/lasso_coef_RR.xlsx")

# lasso models with cross validation for BdgT -----
lasso_coef_bdgt <- c("Intercept",names(dplyr::select(lakes73clean,!c("BdgT")))) %>% as.data.frame() %>% setNames("param")

for(z in 1:M){
  lasso_coef_bdgt %>% tibble::add_column(z = NA)
  data <- complete(RRmice,z)
  coef_fit <- cv.glmnet(formula = BdgT~.,data=data) %>% coef(s="lambda.min") %>% summary()
  for (w in 1:dim(coef_fit)[1]){
    lasso_coef_bdgt[coef_fit[w,1],z+1] <- coef_fit[w,3]
  }
}

lasso_coef_bdgt$pooled <- NA
lasso_coef_bdgt$nb <- NA
for (x in 1:length(lasso_coef_bdgt$param)){
  estimates <- lasso_coef_bdgt[x,] %>% dplyr::select(c(2:M+1)) %>% as.numeric()
  lasso_coef_bdgt$pooled[x] <- mean(estimates,na.rm = T)
  lasso_coef_bdgt$nb[x] <- M - sum(is.na(estimates))
}

write_xlsx(lasso_coef_bdgt,"8.version_control/lasso_coef_BdgT.xlsx")


ggplot(lakes73clean)+geom_point(aes(x=BdgT,y=2.73+2.82*DN+0.15*pH))
lm(BdgT~DN+pH,data=lakes73clean) %>% summary()
