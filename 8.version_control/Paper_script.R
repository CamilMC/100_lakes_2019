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

source("5.100_lakes/100lakes_analysis_functions.R")
norgemap <- st_read("3.Maps/TM_WORLD_BORDERS-0.3.shp")

#load the world map

library("rnaturalearth")
library(rgeos)
world <- ne_countries(scale = "medium", returnclass = "sf")

# load data -----


lakes73 <- read.csv("8.version_control/rr100lakes.csv")

useless.params <- c("X","Survey", "NVE_number","SR")
source("8.version_control/bacteria.R")


# separation of eutrophic/oligotrophic/dystrophic lakes
lakes73$DNOM <- NA
for(i in 1:length(lakes73$DOC)){
  if (lakes73$DOC[i] > mean(lakes73$DOC)){
    lakes73$DNOM[i] <- 1
    } else {
    lakes73$DNOM[i] <- 0
    }
}

lakes73$nuts <- NA
for(i in 1:length(lakes73$DN)){
  if ((lakes73$DN[i]+lakes73$DP[i]) > (mean(lakes73$DN+lakes73$DP))){
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
  scale_color_manual(values=colors_cluster)+
  theme(legend.position = c(.8,.15))+
  xlim(5,15)+ylim(58,62)

ggplot(lakes73)+geom_point(aes(x=DOC,y=RR,col=trophic))
ggplot(lakes73)+geom_point(aes(x=SUVA,y=RR,col=trophic))
ggplot(lakes73)+geom_point(aes(x=SR,y=RR,col=trophic))
ggplot(lakes73)+geom_point(aes(x=DN,y=RR,col=trophic))
ggplot(lakes73)+geom_point(aes(x=DP,y=RR,col=trophic))
ggplot(lakes73)+geom_point(aes(x=DN+DP,y=RR,col=trophic))
ggplot(lakes73)+geom_point(aes(x=DN+DP,y=DOC,col=trophic))
ggplot(lakes73)+geom_point(aes(x=CP,y=RR,col=trophic))
ggplot(lakes73)+geom_point(aes(x=CN,y=RR,col=trophic))



ggplot()+geom_sf(data=norgemap)+
  geom_point(data=lakes73,aes(x=Long,y=Lat,col=trophic),size=4)+
  geom_point(data=filter(lakes73,RR > 10),aes(x=Long,y=Lat,col=trophic),size = 6)+
  theme_void(base_size = 26)+labs(x="",y="",col="trophic")+
  scale_color_manual(values=colors_cluster)+
  theme(legend.position = c(.75,.15))+
  xlim(5,15)+ylim(58,63)

# separation of eutrophic/oligotrophic/dystrophic lakes based on Nünberg (1996)
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

# dataframes with less parameters

lakes73num <- lakes73 %>% select(!c("X","Lake_ID","DNAsampleID","NIVA_station_ID","NVE_number","Lake_name","CBA_sample_date","NIVA_sample_date","SARvis","Survey","incub_date","trophic","DNOM","nuts")) %>% sapply(as.numeric) %>% as.data.frame()


# pca

lakes73num$trophic <- NULL
lakes73pca <- lakes73num %>% fill.na() %>% select(c("DN","DP","DOC","CN","CP","SUVA","SR","Ca","Fe","H","RR","OD","BdgT","O2","CO2","cells_counts.mL"))
rownames(lakes73pca) <- lakes73$Lake_ID
pca73lakes <- prcomp(~.,data=lakes73pca,scale.=T,center=T)
pcafviz <- fviz_pca_biplot(pca73lakes, repel=T,label = c("var"),col.var = "contrib",col.ind = "contrib",
                           gradient.cols =  viridis(3,end=0.8,direction = -1),arrowsize=1,labelsize=6,title="")+ theme_light(base_size=24) 
pcafviz23 <- fviz_pca_biplot(pca73lakes, axes = c(2,3),repel=T,label = c("var"),col.var = "contrib",col.ind = "contrib",
                             gradient.cols =  viridis(3,end=0.8,direction = -1),arrowsize=1,labelsize=6,title="")+ theme_light(base_size=24) 
plot(pcafviz)
plot(pcafviz23)

#without 12777
lakes72 <- filter(lakes73,DOC<40)
lakes72pca <- lakes72 %>% fill.na() %>% select(c("DN","DP","DOC","CN","CP","SUVA","SR","Ca","Fe","H","RR","OD","BdgT","O2","CO2","cells_counts.mL")) %>% sapply(as.numeric) %>% as.data.frame()

rownames(lakes72pca) <- lakes72$Lake_ID
pca72lakes <- prcomp(~.,data=lakes72pca,scale.=T,center=T)
pcafviz <- fviz_pca_biplot(pca72lakes, repel=T,label = c("var"),col.var = "contrib",col.ind = "contrib",
                           gradient.cols =  viridis(3,end=0.8,direction = -1),arrowsize=1,labelsize=6,title="")+ theme_light(base_size=24) 
pcafviz23 <- fviz_pca_biplot(pca72lakes, axes = c(2,3),repel=T,label = c("var"),col.var = "contrib",col.ind = "contrib",
                             gradient.cols =  viridis(3,end=0.8,direction = -1),arrowsize=1,labelsize=6,title="")+ theme_light(base_size=24) 
plot(pcafviz)
plot(pcafviz23)


# pca of outliers
lakes8pca <- lakes73num %>% fill.na() %>% filter(RR > 10) %>% select(c("CN","CP","SUVA","Ca","Fe","H","RR","OD","BdgT","O2","CO2","cells_counts.mL"))

pca8lakes <- prcomp(~.,data=lakes8pca,scale.=T,center=T)
pcafviz <- fviz_pca_biplot(pca8lakes, repel=T,label = "var",col.var = "contrib",col.ind = "contrib",
                           select.var = list(contrib = 15), 
                           gradient.cols =  viridis(3,end=0.8,direction = -1),arrowsize=1,labelsize=6,title="")+ theme_light(base_size=24) 
plot(pcafviz)

pcafviz23 <- fviz_pca_biplot(pca8lakes, axes = c(2,3),repel=T,label = "var",col.var = "contrib",col.ind = "contrib",
                             select.var = list(contrib = 15),
                             gradient.cols =  viridis(3,end=0.8,direction = -1),arrowsize=1,labelsize=6,title="")+ theme_light(base_size=24) 
plot(pcafviz23)

# corrplot

corlakes73 <- cor(lakes73num, use="pairwise.complete.obs",method = c("pearson"))
write_xlsx(as.data.frame(corlakes73),"8.version_control/corlakes73.xlsx")
corlakes73_pvalues <- corrplot::cor.mtest(lakes73num,conf.level=0.95)
corrplot::corrplot(corlakes73,method = c("number"),type=c("upper"),tl.cex = 0.6,number.cex=0.45,p.mat = corlakes73_pvalues$p,insig = "blank")

print.cor.signif2(lakes73,"RR")

# boxplot 
lakes73num$trophic <- lakes73$trophic

pdf("8.version_control/boxplots_lakes73.pdf")

for (i in names(lakes73num)){
  print(ggplot(lakes73num)+geom_boxplot(aes(x=trophic,y=lakes73num[,i]),outlier.shape = NA)+
          geom_jitter(aes(x=trophic,y=lakes73num[,i],col=trophic),size=2)+
          ylab(i)+xlab("")+
          scale_color_viridis_d(end = 0.8)+
          theme_light(base_size=20))
}


dev.off()

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

# clusters -----

csel <- c("Lake_ID","DOC","TN","TP","a410")
lakes73_clusters <- lakes73[csel] 
lakes73_clusters <- apply(lakes73_clusters,2,as.numeric) %>% as.data.frame() 

clust70c <- lakes73_clusters[,-which(names(lakes73_clusters) %in% c("Lake_ID","RR.DOC","RR","RR.OD","OD","BdgT"))]
rownames(clust70c) <- lakes73_clusters$Lake_ID

dd <- dist(scale(clust70c), method = "euclidian")
hc <- hclust(dd, method = "ward.D2")

dend70c <- as.dendrogram(hc)

# optimal number of clusters -----
lakes73_scaled <- scale(clust70c) %>% as.data.frame()
lakes73_scaled <- fill.na(lakes73_scaled)

# elbow method
fviz_nbclust(lakes73_scaled, kmeans, method = "wss") 

#silhouette method
fviz_nbclust(lakes73_scaled, kmeans, method = "silhouette")


fviz_nbclust(lakes73_scaled, kmeans, method = "gap_stat")


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


# correLations ----
pdf("5.100_lakes/all.parameters.correLations.pdf",width=15,height = 8)
sapply(names(lakes73),print.cor.signif2, df=lakes73)
dev.off()

print.cor.signif2(lakes73,"RR.OD")
print.cor.signif2(lakes73,"OD.DOC")
print.cor(lakes73,"RR")
print.cor.signif2(lakes73,"OD")
print.cor.signif2(lakes73,"RR.OD")
print.cor.signif2(lakes73,"DOC")
print.cor.signif2(lakes73,"BdgT")
print.cor.signif2(lakes73,"SAR")
print.cor.signif2(lakes73,"SUVA")

lakes73num$bdgs <- lakes73num$RR / lakes73num$cells_counts.mL
print.cor.signif(lakes73num,"bdgs")

lakes73num$RRsc <- scale(lakes73num$RR,center = T)
lakes73$RRsc <- scale(lakes73$RR,center = T)
print.cor.signif(lakes73num,"RRsc")

ggplot(lakes73)+geom_point(aes(x=log(TOC),y=RRsc,col=trophic))
ggplot(lakes73)+geom_point(aes(x=log(DOC),y=RRsc,col=trophic))
ggplot(lakes73)+geom_point(aes(x=log(DN),y=RRsc,col=trophic))

lakes73num$logRR <- log(lakes73num$RR)
lakes73$logRR <- log(lakes73$RR)

lakes73num$logP <- log(lakes73num$DP)
lakes73num$logN <- log(lakes73num$DN)
log73num$logC <- log(lakes73num$DOC)

print.cor.signif(lakes73num,"logRR")
print.cor.signif(lakes73num,"bdgs")
print.cor.signif(lakes73num,"cells_counts.mL")

ggplot(lakes73)+geom_point(aes(x=cells_counts.mL, y=a600))
ggplot(lakes73)+geom_point(aes(x=logRR, y=a600,col=trophic))



# correLataions - figures ----

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






# boxplots -----
pdf("5. 100 lakes/boxplots_clusters_2.pdf")

for (i in names(lakes73)){
  print(ggplot(lakes73)+geom_boxplot(aes(x=cluster,y=lakes73[,i],group=cluster))+ylab(i)+
          geom_text(aes(x=cluster,y=lakes73[,i], label =Sample_ID),position = "jitter"))
}


dev.off()





# extra plots ----


ggplot(lakes73,aes(x=DOC,y=RR.OD))+geom_point()+geom_smooth(method = "gam",formula = y~s(x))+theme_light(base_size=24)
ggplot(lakes73,aes(x=DN,y=RR.OD))+geom_point()+geom_smooth(method = "gam",formula = y~s(x))+theme_light(base_size=24)
ggplot(lakes73,aes(x=CN,y=RR.OD))+geom_point()+geom_smooth(method = "gam",formula = y~s(x))+theme_light(base_size=24)
ggplot(lakes73,aes(x=pH,y=RR.OD))+geom_point()+geom_smooth(method = "gam",formula = y~s(x))+theme_light(base_size=24)


ggplot(lakes73,aes(x=DOC,y=RR))+geom_point()+geom_smooth(method = "gam",formula = y~s(x))+theme_light(base_size=24)
ggplot(lakes73,aes(x=DN,y=RR))+geom_point()+geom_smooth(method = "gam",formula = y~s(x))+theme_light(base_size=24)
ggplot(lakes73,aes(x=CN,y=RR))+geom_point()+geom_smooth(method = "gam",formula = y~s(x))+theme_light(base_size=24)

