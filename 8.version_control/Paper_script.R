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

source("5. 100 lakes/100lakes_analysis_functions.R")
norgemap <- st_read("3. Maps/TM_WORLD_BORDERS-0.3.shp")


# load data -----



lakes70 <- read.csv("5. 100 lakes/lakes70.csv")

useless.params <- c("X","NVE_number","DNA_ug.mL","filtered_mL","SAR","s_350_400",
                    "Lake_Area","ionic_strength_.mol.L","cluster","Al","CO2_nM.h","N2_nM.h","O2_nM.h","CH4_nM.h","N2O_nM.h",
                    "TCN","TCP","week")
source("8. version control/bacteria.R")


# maps -----

pdf("5. 100 lakes/Vmax_maps.pdf",width=10,height=10)

ggplot()+geom_sf(data=norgemap)+
  geom_point(data=lakes70,aes(x=long,y=lat,col=Vmax),size=5)+
  geom_point(data=filter(lakes70,Vmax >8),aes(x=long,y=lat,size=Vmax))+scale_size(limits=c(8,30),breaks=c(10,20,30),range=c(4,6))+
  theme_void(base_size = 24)+labs(x="",y="",col=expression(atop("Vmax",paste("(",mu,"mol/L/h)"))),size="outliers")+
  scale_color_gradientn(colors = viridis(9,direction = -1,end=0.9),limits=c(0,8))+
  #geom_text(data=filter(lakes70,Vmax>8),aes(x=long,y=lat,label=round(Vmax,0)),nudge_y=0.1)+
  guides(color=guide_legend(order=1),size = guide_legend(order=2))+
  theme(legend.position = c(.87,.45))+
  xlim(5,15)+ylim(58,63)

ggplot()+geom_sf(data=norgemap)+
  geom_point(data=lakes70,aes(x=long,y=lat,col=Vmax.LOC),size=5)+
  theme_void(base_size = 24)+labs(x="",y="",col=expression(atop("Vmax/LOC",(h^{-1}))))+
  scale_color_viridis_c(direction = -1, end = 0.9)+
  theme(legend.position = c(.9,.45))+
  xlim(5,15)+ylim(58,63)

ggplot()+geom_sf(data=norgemap)+
  geom_point(data=lakes70,aes(x=long,y=lat,col=Vmax.DOC),size=5)+
  geom_point(data=filter(lakes70,Vmax.DOC >2),aes(x=long,y=lat,size=Vmax.DOC))+scale_size(limits=c(2,10),breaks=c(3,5,10),range=c(4,6))+
  theme_void(base_size = 24)+labs(x="",y="",col=expression(atop("Vmax/DOC",(h^{-1}))),size="outliers")+
  scale_color_gradientn(colors = viridis(6,direction = -1,end=0.9),limits=c(0,2))+
  #geom_text(data=filter(lakes70,Vmax.DOC >2),aes(x=long,y=lat,label=round(Vmax.DOC,2)),nudge_y = -0.09,nudge_x = -0.1)+
  guides(color=guide_legend(order=1),size = guide_legend(order=2))+
  theme(legend.position = c(.87,.45))+
  xlim(5,15)+ylim(58,63)

ggplot()+geom_sf(data=norgemap)+
  geom_point(data=lakes70,aes(x=long,y=lat,col=LOC.DOC),size=4)+
  theme_void(base_size = 24)+labs(x="",y="",col="%DOC\nconsumed",size="outliers")+
  scale_color_gradientn(colors = viridis(6,direction = -1,end=0.9),limits=c(0,4))+
  geom_point(data=filter(lakes70,LOC.DOC >4),aes(x=long,y=lat,size=LOC.DOC))+scale_size(limits=c(4,20),breaks=c(5,10,15),range=c(4,6))+
#  geom_text(data=filter(lakes70,LOC.DOC >3),aes(x=long,y=lat,label=round(LOC.DOC,2)),nudge_y = -0.09,nudge_x = -0.1)+  
  guides(color=guide_legend(order=1),size = guide_legend(order=2))+
  theme(legend.position = c(.87,.5))+
  xlim(5,15)+ylim(58,63)

ggplot()+geom_sf(data=norgemap)+
  geom_point(data=lakes70,aes(x=long,y=lat,col=LOC),size=4)+
  theme_void(base_size = 24)+labs(x="",y="",col="LOC\n(umol/L)")+
  scale_color_gradientn(colors = viridis(6,direction = -1,end=0.9),limits=c(0,30))+
  geom_text(data=filter(lakes70,LOC >30),aes(x=long,y=lat,label=round(LOC,0)),nudge_y = -0.09,nudge_x = -0.1)+  
  theme(legend.position = c(.9,.7))+
  xlim(5,15)+ylim(58,63)

ggplot()+geom_sf(data=norgemap)+
  geom_point(data=lakes70,aes(x=long,y=lat,col=DN),size=4)+
  geom_point(data=filter(lakes70,DN > 0.6),aes(x=long,y=lat,size= DN))+scale_size(limits=c(0.6,0.9))+
  theme_void(base_size = 24)+labs(x="",y="",col="DN\nmg/L",size="Outliers")+
  scale_color_gradientn(colors = viridis(6,direction = -1,end=0.9),limits=c(0,0.6))+
#  geom_text(data=filter(lakes70,DN > 0.6),aes(x=long,y=lat,label=round(DN,1)),nudge_y = -0.09,nudge_x = -0.1)+  
  theme(legend.position = c(.9,.6))+
  xlim(5,15)+ylim(58,63)

ggplot()+geom_sf(data=norgemap)+
  geom_point(data=lakes70,aes(x=long,y=lat,col=DP),size=4)+
  geom_point(data=filter(lakes70,DP > 10),aes(x=long,y=lat,size= DP))+scale_size(limits=c(10,16))+
  theme_void(base_size = 24)+labs(x="",y="",col="DP\nug/L",size="Outliers")+
  scale_color_gradientn(colors = viridis(6,direction = -1,end=0.9),limits=c(0,10))+
  #  geom_text(data=filter(lakes70,DN > 0.6),aes(x=long,y=lat,label=round(DN,1)),nudge_y = -0.09,nudge_x = -0.1)+  
  theme(legend.position = c(.9,.6))+
  xlim(5,15)+ylim(58,63)

ggplot()+geom_sf(data=norgemap)+
  geom_point(data=lakes70,aes(x=long,y=lat,col=CN),size=4)+
  geom_point(data=filter(lakes70,CN > 50),aes(x=long,y=lat,size= CN))+scale_size(limits=c(50,150),range=c(4,6),breaks=c(50,100,150))+
  theme_void(base_size = 24)+labs(x="",y="",col="CN",size="Outliers")+
 # scale_color_gradientn(colors = viridis(6,direction = -1,end=0.9))+
  scale_color_gradientn(colors = viridis(6,direction = -1,end=0.9),limits=c(8,50))+
 # geom_text(data=filter(lakes70,DN > 0.6),aes(x=long,y=lat,label=round(DN,1)),nudge_y = -0.09,nudge_x = -0.1)+ 
  guides(color=guide_legend(order=1),size = guide_legend(order=2))+
  theme(legend.position = c(.9,.6))+
  xlim(5,15)+ylim(58,63)

ggplot()+geom_sf(data=norgemap)+
  geom_point(data=lakes70,aes(x=long,y=lat,col=DOC),size=5)+
  #geom_point(data=filter(lakes70,DOC > 20),aes(x=long,y=lat,size= DOC))+scale_size(limits=c(20,21))+
  theme_void(base_size = 26)+labs(x="",y="",col="DOC\nmg/L",size="Outliers")+
  # scale_color_gradientn(colors = viridis(6,direction = -1,end=0.9))+
  scale_color_gradientn(colors = viridis(6,direction = -1,end=0.9),limits=c(1,20))+
#  geom_text(data=filter(lakes70,DOC > 15),aes(x=long,y=lat,label=round(DOC,1)),nudge_y = -0.1,nudge_x = -0.0)+ 
  guides(color=guide_legend(order=1),size = guide_legend(order=2))+
  theme(legend.position = c(.75,.2))+
  xlim(5,15)+ylim(58,63)

dev.off()


# clusters -----

csel <- c("Sample_ID",
          "temperature","pH","conductivity","alkalinity",
          "DOC","DN","DP","TOC","TN","TP",
          "CN","CP",
          "CO2","CH4","O2",
          "cells",
          "Cl","SO4",
          "SUVA","SR",
          "Altitude","Mg","Ca_NIVA","Na_NIVA","K","Fe","Al_NIVA",
          "Ni","Cu","As","V","Co")
lakes70_clusters <- lakes70[csel] 
clust70c <- lakes70_clusters[,-which(names(lakes70_clusters) %in% c("X","Sample_ID","long","lat","cluster"))]
rownames(clust70c) <- lakes70_clusters$Sample_ID

dd <- dist(scale(clust70c), method = "euclidian")
hc <- hclust(dd, method = "ward.D2")

dend70c <- as.dendrogram(hc)

# optimal number of clusters -----
lakes70_scaled <- scale(lakes70_clusters) %>% as.data.frame()
lakes70_scaled <- fill.na(lakes70_scaled)

# elbow method
fviz_nbclust(lakes70_scaled, kmeans, method = "wss") 

#silhouette method
fviz_nbclust(lakes70_scaled, kmeans, method = "silhouette")


fviz_nbclust(lakes70_scaled, kmeans, method = "gap_stat")


# creates clusters ----
nb_clusters <- 4
lakes70$cluster <- cutree(dend70c,nb_clusters)
colors_cluster <- viridis(nb_clusters,end = 0.8,direction = -1)


dendplotc <- dend70c %>% set("labels_col",value = colors_cluster,k=nb_clusters) %>%
  set("labels_cex",0.8) %>%
  set("branches_k_color", value = colors_cluster,k=nb_clusters)
plot(dendplotc,horiz = T)

png("5. 100 lakes/dendrogram.png",height = 800, width=300)
plot(dendplotc,horiz = T)
dev.off()


# graphs and map of clusters -----
pdf("5. 100 lakes/bdg_by_cluster.pdf",width=10,height=8)
ggplot(lakes70,aes(y=Vmax.LOC,x=cluster))+geom_boxplot(aes(group=cluster),outlier.shape = NA)+geom_jitter(aes(col=DOC),width=0.3,size=3)+
  scale_color_viridis_c(direction = -1,end=0.9)+theme_light(base_size=26)+ylab("Vmax/LOC")
ggplot(lakes70,aes(y=Vmax,x=cluster))+geom_boxplot(aes(group=cluster),outlier.shape = NA)+geom_jitter(aes(col=DOC),width=0.3,size=3)+
  scale_color_viridis_c(direction = -1,end=0.9)+theme_light(base_size=26)+ylab("Vmax")
ggplot(lakes70,aes(y=Vmax.DOC,x=cluster))+geom_boxplot(aes(group=cluster),outlier.shape = NA)+geom_jitter(aes(col=DOC),width=0.3,size=3)+
  scale_color_viridis_c(direction = -1,end=0.9)+theme_light(base_size=26)+ylab("Vmax/DOC")

dev.off()

ggplot()+geom_sf(data=norgemap)+
  geom_point(data=lakes70,aes(x=long,y=lat,col=as.character(cluster)),size=6)+
  theme_void(base_size = 26)+labs(x="",y="",col="cluster")+
  scale_color_manual(values=colors_cluster)+
  theme(legend.position = c(.75,.15))+
  xlim(5,15)+ylim(58,63)
ggsave("5. 100 lakes/map_clusters.png",width=10,height = 10)

hist(lakes70$log.Vmax)

# ANOVA -----

a1 <- aov(Vmax.LOC~as.factor(cluster),data=lakes70_clusters)
posthoc <- TukeyHSD(x=a1,conf.level=0.95)

aov_summary <- as.data.frame(names(lakes70_clusters))
aov_summary$pvalue <- NA
aov_summary$effect_size <- NA
aov_summary$ptt12 <- NA
aov_summary$ptt13 <- NA
aov_summary$ptt23 <- NA

for (i in 1:length(names(lakes70_clusters))){
  param <- names(lakes70_clusters)[i]
  aovt <- aov(lakes70_clusters[,param]~lakes70_clusters$cluster)
  pval <- summary(aovt)[[1]]$`Pr(>F)`[1]
  es <- eta_squared(aovt)$Eta_Sq_partial
  ptt <- pairwise.t.test(lakes70_clusters[,param],lakes70_clusters$cluster,p.adjust.method = "bonferroni")
  aov_summary$pvalue[i] <- pval
  aov_summary$effect_size[i] <- es
  aov_summary$ptt12[i] <- ptt$p.value[1]
  aov_summary$ptt13[i] <- ptt$p.value[2]
  aov_summary$ptt23[i] <- ptt$p.value[4]
}

names(aov_summary)[1] <- "param"
write_xlsx(aov_summary,path="5. 100 lakes/anovas_clusters_fdr.xlsx")


ptt <- aov_summary %>% select(c("param","ptt12","ptt13","ptt23","pvalue"))

ptt$pvalue[which(ptt$pvalue > 0.05)] <- NA
ptt$pvalue[which(ptt$pvalue < 0.05)] <- "significant"
ptt$pvalue[which(is.na(ptt$pvalue))] <- "non-significant"

ptt$ptt12[which(ptt$ptt12 > 0.05)] <- NA
ptt$ptt12[which(ptt$ptt12 < 0.05)] <- "Cluster 1 ??? Cluster 2"

ptt$ptt13[which(ptt$ptt13 > 0.05)] <- NA
ptt$ptt13[which(ptt$ptt13 < 0.05)] <- "Cluster 1 ??? Cluster 3"

ptt$ptt23[which(ptt$ptt23 > 0.05)] <- NA
ptt$ptt23[which(ptt$ptt23 < 0.05)] <- "Cluster 2 ??? Cluster 3"


ptt.graph <- ptt %>% filter(!param %in% c("Sample_ID","cluster"))

#ptt.graph <- ptt %>% filter(!param %in% c("V51","Sample_ID","lag_bdg","cluster","Cu63","Ni58","V51","As75"))

ptt.graph$param[which(ptt.graph$param == "Vmax.LOC")] <- "Biodeg_LOC"
ptt.graph$param[which(ptt.graph$param == "Vmax")] <- "Biodeg"
ptt.graph$param[which(ptt.graph$param == "Vmax.DOC")] <- "Biodeg_s"
ptt.graph$param[which(ptt.graph$param == "LOC.DOC")] <- "%DOM"
ptt.graph$param[which(ptt.graph$param == "pH_lab")] <- "pH"
ptt.graph$param[which(ptt.graph$param == "alk_meq.L")] <- "Alkalinity"
ptt.graph$param[which(ptt.graph$param == "cond_µS.m")] <- "Conductivity"
ptt.graph$param[which(ptt.graph$param == "CO2_µM")] <- "CO2"
ptt.graph$param[which(ptt.graph$param == "Chloride_mg.L")] <- "Cl"
ptt.graph$param[which(ptt.graph$param == "Sulfate_mg.L")] <- "SO4"
ptt.graph$param[which(ptt.graph$param == "T")] <- "Bdg period"
ptt.graph$param[which(ptt.graph$param == "temp")] <- "Temperature"
ptt.graph$param[which(ptt.graph$param == "CH4_nM")] <- "CH4"
ptt.graph$param[which(ptt.graph$param == "cond_µS.m")] <- "Conductivity"
ptt.graph$param[which(ptt.graph$param == "Al_R")] <- "Al"
ptt.graph$param[which(ptt.graph$param == "O2_µM")] <- "O2"
ptt.graph$param[which(ptt.graph$param == "Ca.y")] <- "Ca"
ptt.graph$param[which(ptt.graph$param == "Na.y")] <- "Na"


x <- c("Biodeg_LOC", "CN","pH","Al",
       "TOC","DOC","Fe","Altitude","Co59","V51","Temperature","Cu63",
       "Na","DN","Mg","K","Ni58",
       "Alkalinity","Conductivity","CO2","TN","TP","Cl-","SO4","Bdg period",
       "SUVA", "As75",
       "CH4","cells_counts.mL","SR","CP",
       "DP","Ca","O2","tmax","Biodeg","Biodeg_s","%DOM",
       "log.biodeg","log.biodegs")
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

k1 <- kruskal.test(DOC~cluster,data=lakes70)
pwt <- pairwise.wilcox.test(lakes70$Vmax,lakes70$cluster,p.adjust.method = "holm")

kw_summary <- as.data.frame(names(lakes70))
kw_summary$pvalue <- NA
#aov_summary$effect_size <- NA
kw_summary$pwt12 <- NA
kw_summary$pwt13 <- NA
kw_summary$pwt23 <- NA

for (i in 1:length(names(lakes70))){
  param <- names(lakes70)[i]
  kwt <- kruskal.test(lakes70[,param]~lakes70$cluster)
  pval <- kwt$p.value
  #es <- eta_squared(aovt)$Eta_Sq_partial
  pwt <- pairwise.wilcox.test(lakes70[,param],lakes70$cluster,p.adjust.method = "holm",paired=FALSE)
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
       "Vmax.LOC", "CN","CP","pH","SR","BdgT",
       "TOC","DOC","Vmax.DOC","LOC.DOC","Altitude",
       "Na","V","As","Ni",
       "DN","Mg","K",
       "Alkalinity","Conductivity","CO2","TN","TP","SO4","DP",
      "Ca","cells",
      "O2","tmax","Vmax",
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
A <- ggplot(lakes70)+geom_boxplot(aes(x=Vmax,y="Biodeg\n(uM/h)"))+theme_light(base_size=15)+labs(x="",y="")
A2 <- ggplot(lakes70)+geom_boxplot(aes(x=LOC,y="LOC\n(uM)"))+theme_light(base_size=15)+labs(x="",y="")
B <- ggplot(lakes70)+geom_boxplot(aes(x=Vmax.DOC,y="Biodeg_s\n(h-1)"))+theme_light(base_size=15)+labs(x="",y="")
C <- ggplot(lakes70)+geom_boxplot(aes(x=Vmax.LOC,y="Biodeg_LOC\n(h-1)"))+theme_light(base_size=15)+labs(x="",y="")
gA <- ggplotGrob(A)
gA2 <- ggplotGrob(A2)
gB <- ggplotGrob(B)
gC <- ggplotGrob(C)
grid::grid.newpage()
grid::grid.draw(cbind(rbind(gA,gA2),rbind(gB,gC)))

D <- ggplot(lakes70)+geom_boxplot(aes(x=DOC,y="DOC (mg/L)"))+theme_light(base_size=15)+labs(x="",y="")
E <- ggplot(lakes70)+geom_boxplot(aes(x=DN,y="DN (mg/L)"))+theme_light(base_size=15)+labs(x="",y="")
F <- ggplot(lakes70)+geom_boxplot(aes(x=DP,y="DP (ug/L)"))+theme_light(base_size=15)+labs(x="",y="")
G <- ggplot(lakes70)+geom_boxplot(aes(x=TOC,y="TOC (mg/L)"))+theme_light(base_size=15)+labs(x="",y="")
H <- ggplot(lakes70)+geom_boxplot(aes(x=TN,y="TN (mg/L)"))+theme_light(base_size=15)+labs(x="",y="")
I <- ggplot(lakes70)+geom_boxplot(aes(x=TP,y="TP (ug/L)"))+theme_light(base_size=15)+labs(x="",y="")
J <- ggplot(lakes70)+geom_boxplot(aes(x=DOC/DN,y="C:N (dissolved)"))+theme_light(base_size=15)+labs(x="",y="")
K <- ggplot(lakes70)+geom_boxplot(aes(x=TOC/TN,y="C:N (total)"))+theme_light(base_size=15)+labs(x="",y="")
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

L <- ggplot(lakes70)+geom_boxplot(aes(x=pH_lab,y="pH"))+theme_light(base_size=15)+labs(x="",y="") 
gL <- ggplotGrob(L)
M <- ggplot(lakes70)+geom_boxplot(aes(x=alk_meq.L,y="alkalinity (meq/L)"))+theme_light(base_size=15)+labs(x="",y="")
gM <- ggplotGrob(M)
N <- ggplot(lakes70)+geom_boxplot(aes(x=Altitude,y="Altitude (m)"))+theme_light(base_size=15)+labs(x="",y="")
grid::grid.newpage()
grid::grid.draw(rbind(gL,gM))

N <- ggplot(lakes70)+geom_boxplot(aes(x=Altitude,y="Altitude (m)"))+theme_light(base_size=15)+labs(x="",y="")
gN <- ggplotGrob(N)

grid::grid.newpage()
grid::grid.draw(cbind(rbind(gD,gE,gJ),rbind(gF,gL,gN)))


ggplot(lakes70)+geom_point(aes(x=LOC,y=Vmax,col=DOC))+scale_color_viridis_c(direction = -1, end = 0.9)+theme_light(base_size=20)

dev.off()

summary_lakes70 <- do.call(cbind,lapply(lakes70,summary)) %>% as.data.frame()
summary_lakes70[,1] <- rownames(summary_lakes70)
write_xlsx(summary_lakes70,"5. 100 lakes/summary_lakes70.xlsx")



# Charge balance with ions -----
lakes70_ions <- lakes70 %>% select(c("Sample_ID","Ca.y","Mg","Na.y","K","Al_R",
           "Sulfate_mg.L","Nitrate_mg.L","Nitrite_mg.L","Chloride_mg.L",
           "alk_meq.L","cond_µS.m","CO2_µM","pH_lab","TOC"))
lakes70_ions$Ca_meq.L <- lakes70_ions$Ca.y/40.078*2 
lakes70_ions$Mg_meq.L <- lakes70_ions$Mg/24.305*2
lakes70_ions$Na_meq.L <- lakes70_ions$Na.y/23.0
lakes70_ions$K_meq.L <- lakes70_ions$K/39.0983
lakes70_ions$Al_meq.L <- lakes70_ions$Al_R*10^(-3)/26.982*3
lakes70_ions$Sulfate_meq.L <- lakes70_ions$Sulfate_mg.L/(32.06+4*15.999)*(-2)
lakes70_ions$Nitrate_meq.L <- lakes70_ions$Nitrate_mg.L/(14+3*15.999)*(-1)
lakes70_ions$Nitrite_meq.L <- lakes70_ions$Nitrite_mg.L/(14+2*15.999)*(-1)
lakes70_ions$Chloride_meq.L <- lakes70_ions$Chloride_mg.L/(35.45)*(-1)
lakes70_ions$Phosphate_meq.L <- lakes70_ions$Phosphate_mg.L/(30.974+4*15.999)*(-3)

# carbonate speciation
Kh <- 0.0316227766016838
K1 <- 5.01187233627272E-07
K2 <- 5.62341325190349E-11

lakes70_ions$co2_mol.L <- lakes70_ions$CO2_µM*10^(-6)
lakes70_ions$h2co3 <- lakes70_ions$co2_mol.L* Kh
lakes70_ions$hco3 <- lakes70_ions$h2co3 * K1 / 10^(-lakes70_ions$pH_lab)
lakes70_ions$co3 <- lakes70_ions$hco3 * K2 / 10^(-lakes70_ions$pH_lab)

lakes70_ions$hco3_meq.L <- lakes70_ions$hco3*10^3*(-1)

# organic anions
lakes70_ions$pK_organions <- 0.96 + 0.90 * lakes70_ions$pH_lab - 0.039 * lakes70_ions$pH_lab ^2
lakes70_ions$organions_meq.L <- - (lakes70_ions$TOC*10) * 10^(-lakes70_ions$pK_organions) / (10^(-lakes70_ions$pH_lab)+10^(-lakes70_ions$pK_organions)) / 1000

# balance
lakes70_ions <- lakes70_ions %>% rowwise() %>% 
  mutate(sum_cations = sum(Ca_meq.L,Mg_meq.L,Na_meq.L,K_meq.L,Al_meq.L, na.rm = TRUE))
lakes70_ions <-  lakes70_ions %>% rowwise() %>%
  mutate(sum_anions = sum(Sulfate_meq.L,Nitrate_meq.L,Chloride_meq.L,Nitrite_meq.L,
                          hco3_meq.L,organions_meq.L,na.rm = TRUE)) 
lakes70_ions$balance <- lakes70_ions$sum_cations + lakes70_ions$sum_anions
lakes70_ions$EB <- (lakes70_ions$sum_cations + lakes70_ions$sum_anions)/
  (lakes70_ions$sum_cations - lakes70_ions$sum_anions) * 100

#bar chart
lakes70_cations_melt <- lakes70_ions %>% select(c("Sample_ID", "Ca_meq.L","Mg_meq.L","Na_meq.L","K_meq.L")) %>%
  melt(id.vars = "Sample_ID",measure.vars = c("Ca_meq.L","Mg_meq.L","Na_meq.L","K_meq.L"),
       variable.name = "ion", value.name = "concentration_meq.L")
lakes70_cations_melt$sample_axis <- paste(lakes70_cations_melt$Sample_ID,".1",sep="")
lakes70_cations_melt$type <- "cation"

lakes70_anions_melt <- lakes70_ions %>% 
  select(c("Sample_ID", "Sulfate_meq.L","Nitrate_meq.L","Chloride_meq.L","hco3_meq.L","organions_meq.L")) %>%
  melt(id.vars = "Sample_ID",measure.vars = c("Sulfate_meq.L","Nitrate_meq.L","Chloride_meq.L","hco3_meq.L","organions_meq.L"),
       variable.name = "ion", value.name = "concentration_meq.L")
lakes70_anions_melt$concentration_meq.L <- lakes70_anions_melt$concentration_meq.L * (-1)
lakes70_anions_melt$sample_axis <- paste(lakes70_anions_melt$Sample_ID,".2",sep="")
lakes70_anions_melt$type <- "anion"


lakes70_ions_melt <- rbind(lakes70_cations_melt,lakes70_anions_melt)
lakes70_ions_melt$type <- factor(lakes70_ions_melt$type,levels = c("cation","anion"))

ggplot(lakes70_ions_melt,aes(x=type,y=concentration_meq.L,fill=ion))+geom_bar(stat = "identity")+
  scale_fill_viridis_d(option = "viridis",begin = 0.2)+ facet_wrap(Sample_ID~.,ncol = 35,strip.position = "bottom")+
  theme_minimal(base_size = 15)+
  theme(axis.text.x = element_blank(),
        legend.title = element_blank(),panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        panel.grid.major.y = element_line(colour="grey95"),
        panel.grid.minor.y = element_line(colour="grey95"))+
  labs(x="Sample",y="meq/L")
ggsave("5. 100 lakes/charge_balance_ions_niva.png",width = 20, height = 10)

# correlations ----
print.cor.signif2(lakes70,"Vmax.LOC")
print.cor.signif2(lakes70,"LOC.DOC")
print.cor(lakes70,"Vmax")
print.cor.signif2(lakes70,"LOC")
print.cor.signif2(lakes70,"Vmax.DOC")
print.cor.signif2(lakes70,"DOC")

# correlataions - figures ----

pdf("5. 100 lakes/Vmax_corr.pdf",height = 10,width=15)
cor_param_Vmax <- c("Sample_ID","LOC","Vmax","CN","CP") 
cor_Vmax <- lakes70[cor_param_Vmax] %>% setNames(c("Sample_ID","LOC","Vmax","C:N","C:P"))
print.cor.signif2(cor_Vmax,"Vmax",title=F,y1=-0.5)

cor_param_Vmax.DOC <- c("Sample_ID","LOC","Vmax","s_275_295","SR","DN","Al_NIVA","SUVA","DOC","Vmax.DOC","LOC.DOC")
cor_Vmax.DOC <- lakes70[cor_param_Vmax.DOC] %>% setNames(c("Sample_ID","LOC","Vmax","s_275_295","SR","DN","Al","SUVA","DOC","Vmax/DOC","%DOM"))
print.cor.signif2(cor_Vmax.DOC,"Vmax/DOC",title=F,y1=-0.5)
print.cor.signif2(cor_Vmax.DOC,"%DOM",title=F,y1=-0.5,y2=1.3)

cor_param_Vmax.LOC <- c("Sample_ID","Al_NIVA","Mg","K","DN","alkalinity","pH","Vmax.LOC","CN")
cor_Vmax.LOC <- lakes70[cor_param_Vmax.LOC] %>% setNames(c("Sample_ID","Al","Mg","K","DN","alkalinity","pH","Vmax/LOC","C:N"))
print.cor.signif2(cor_Vmax.LOC,"Vmax/LOC",title = F,y1=-0.7,y2=0.5)

cor_param_DOC <- c("DOC","TOC","V","Al","SUVA","Cu","Na","Fe","Co","Pb","Zn","cells","Br","DN","TN","CP","CN","K","conductivity","Ni","Mg","CH4","Vmax.DOC","LOC.DOC","temperature","SR","pH","Altitude")
cor_DOC <- lakes70[cor_param_DOC]
print.cor.signif2(cor_DOC,"DOC",title = F)

dev.off()

# PCA -----
pca_params <- c("Altitude","pH","SR","Vmax.LOC","DN","SUVA","DOC","alkalinity","Vmax","Vmax.DOC","CN","CP")
lakes70pca <- lakes70[pca_params]
pca70lakes <- prcomp(~.,data=lakes70pca,scale.=T,center=T)
pcafviz <- fviz_pca_biplot(pca70lakes, repel=T,label = "var",col.var = "contrib",col.ind = "contrib",
                           gradient.cols =  viridis(3,end=0.8,direction = -1),arrowsize=1,labelsize=6,title="")+ theme_light(base_size=24) 
pcafviz23 <- fviz_pca_biplot(pca70lakes, axes = c(2,3),repel=T,label = "var",col.var = "contrib",col.ind = "contrib",
                             gradient.cols =  viridis(3,end=0.8,direction = -1),arrowsize=1,labelsize=6,title="")+ theme_light(base_size=24) 
plot(pcafviz)
ggsave("5. 100 lakes/Biodeg_pca.png",height = 8, width = 10)


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

for (i in names(lakes70)){
  print(ggplot(lakes70)+geom_boxplot(aes(x=cluster,y=lakes70[,i],group=cluster))+ylab(i)+
          geom_text(aes(x=cluster,y=lakes70[,i], label =Sample_ID),position = "jitter"))
}


dev.off()

pdf("5. 100 lakes/boxplots_lakes70.pdf")
for (i in names(clust70c)){
  print(ggplot(clust70c)+geom_boxplot(aes(x=i,y=clust70c[,i]))+ylab(i)+xlab("")+
          theme_light(base_size=20))
}

for (i in names(lakes70)){
  print(ggplot(lakes70)+geom_boxplot(aes(x=cluster,y=lakes70[,i],group=cluster))+ylab(i)+xlab("")+
          theme_light(base_size=20))
}

dev.off()



# extra plots ----


ggplot(lakes70,aes(x=DOC,y=Vmax.LOC))+geom_point()+geom_smooth(method = "gam",formula = y~s(x))+theme_light(base_size=24)
ggplot(lakes70,aes(x=DN,y=Vmax.LOC))+geom_point()+geom_smooth(method = "gam",formula = y~s(x))+theme_light(base_size=24)
ggplot(lakes70,aes(x=CN,y=Vmax.LOC))+geom_point()+geom_smooth(method = "gam",formula = y~s(x))+theme_light(base_size=24)
ggplot(lakes70,aes(x=pH,y=Vmax.LOC))+geom_point()+geom_smooth(method = "gam",formula = y~s(x))+theme_light(base_size=24)


ggplot(lakes70,aes(x=DOC,y=Vmax))+geom_point()+geom_smooth(method = "gam",formula = y~s(x))+theme_light(base_size=24)
ggplot(lakes70,aes(x=DN,y=Vmax))+geom_point()+geom_smooth(method = "gam",formula = y~s(x))+theme_light(base_size=24)
ggplot(lakes70,aes(x=CN,y=Vmax))+geom_point()+geom_smooth(method = "gam",formula = y~s(x))+theme_light(base_size=24)

