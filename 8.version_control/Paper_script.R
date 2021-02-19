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
library(miceadds)
library(glmnet)
library(glmnetUtils)
library(rstatix)

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

ggplot(lakes73)+geom_point(aes(x=RR,y=BdgT,col=OD),size = 4,shape = 5)+theme_light(base_size=20)+scale_color_viridis_c(end = 0.8, direction = -1)+
  labs(x="Respiration rate (umol/h)",y="Biodegradation period (h)", col = "Oxygen\ndemand\n(umol)")
ggsave("BdgTvsRR.png",width = 10,height = 6)

boxplot(lakes73$RR)
boxplot(lakes73$BdgT)

g <- ggplot()+geom_boxplot(data=lakes73,aes(x="",y=RR))
outliers <- ggplot_build(g)[["data"]][[1]][["outliers"]] %>% as.data.frame()
outliers_ID <-  filter(lakes73, lakes73$RR %in% outliers[[1]]) %>% dplyr::select(Lake_name)
outliers$Lake_ID <- outliers_ID$Lake_name
names(outliers) <- c("RR","Lake_name")
g <- ggplot()+geom_boxplot(data=lakes73,aes(x="",y=RR),outlier.size = 3, width = 1.5)+theme_light(base_size = 20)+xlab("")+ylab("RR (umol/h)")+
  geom_text_repel(data=outliers,aes(x="",y=outliers[,"RR"],label=Lake_name),position = position_dodge(width=2),size = 6, col = "steelblue4",segment.size = 1)
print(g)
ggsave("boxplot_RR.png",width=4,height=6)

h <- ggplot()+geom_boxplot(data=lakes73,aes(x="",y=BdgT))
outliers <- ggplot_build(h)[["data"]][[1]][["outliers"]] %>% as.data.frame()
outliers_ID <-  filter(lakes73, lakes73$BdgT %in% outliers[[1]]) %>% dplyr::select(Lake_name)
outliers$Lake_ID <- outliers_ID$Lake_name
names(outliers) <- c("BdgT","Lake_name")
h <- ggplot()+geom_boxplot(data=lakes73,aes(x="",y=BdgT),outlier.size = 3,width=1.5)+theme_light(base_size = 20)+xlab("")+ylab("BdgT (h)")+
  geom_text_repel(data=outliers,aes(x="",y=outliers[,"BdgT"],label=Lake_name),position = position_dodge(width=2),size = 6, col = "steelblue4",segment.size = 1)
print(h)
ggsave("boxplot_BdgT.png",width=4,height=6)
       
boxplots <- function(df,coef_out){
  for(i in names(df)[-1]){
    g <- ggplot()+geom_boxplot(data=df,aes(x=i,y=df[,i]),coef=coef_out)
    outliers <- ggplot_build(g)[["data"]][[1]][["outliers"]] %>% as.data.frame()
    outliers_ID <-  filter(df, df[,i] %in% outliers[[1]]) %>% dplyr::select(Sample_ID)
    outliers$Sample_ID <- outliers_ID$Sample_ID
    names(outliers) <- c(i,"Sample_ID")
    if(length(outliers$Sample_ID) == 0){
      g <- ggplot()+geom_boxplot(data=df,aes(x=i,y=df[,i]),outlier.size = 4,coef = coef_out)+theme_light(base_size = 20)+xlab("")+ylab(i)
      print(g)
    } else {
      g <- ggplot()+geom_boxplot(data=df,aes(x=i,y=df[,i]),outlier.size = 4,coef = coef_out)+theme_light(base_size = 20)+xlab("")+ylab(i)+
        geom_text_repel(data=outliers,aes(x=i,y=outliers[,i],label=Sample_ID),position = position_dodge(width=1),size = 6, col = "red",segment.size = 1)
      print(g)
    }
    out[[i]] <- outliers
  }
}


# -----
# distribution of data ----
pdf("8.version_control/distribution.pdf")
hist.df(lakes73num)
dev.off()

# -----
# corrplot -----

corlakes73 <- cor(lakes73num, use="complete.obs",method = c("pearson"))
write_xlsx(as.data.frame(corlakes73),"8.version_control/corlakes73.xlsx")
corlakes73_pvalues <- corrplot::cor.mtest(lakes73num,conf.level=0.95)
corrplot::corrplot(corlakes73,method = c("number"),type=c("upper"),tl.cex = 0.6,number.cex=0.45,p.mat = corlakes73_pvalues$p,insig = "blank")

print.cor.signif2(lakes73,"RR")

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
# maps ------

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
# correLations ----
pdf("5.100_lakes/all.parameters.correLations.pdf",width=15,height = 8)
sapply(names(lakes73),print.cor.signif2, df=lakes73)
dev.off()

# -----
# create multiple imputation datasets to fill NA -----
param2remove <- c("auc","width","OD","Lake_ID","NVE_number","NIVA_date","X","Long","Lat","Altitude","CBA_week","tmax","tmax_h","ox_initial","F","pH_bio","EC_bio","NO2","NO3","DNA","p_N2","p_O2","c_O2","p_CO2","c_CO2","p_CH4","p_N2O","d18O","d2H","Tso","a254","a400","a410","a600","a275","a295","s_275_295","a350","a400","s_350_400","lag_bdg","H","NP","CNP","dist_bact","PO4")
lakes73clean <- lakes73num %>% dplyr::select(!which(names(lakes73num)%in% param2remove))
write_xlsx(lakes73clean,"8.version_control/lakes73clean.xlsx")

M <- 10
set.seed(5)
RRmice <- lakes73clean %>% mice(method = "cart",m = M)
# -----
# cor.test on mice RR -----
predvar <- names(lakes73clean)[-which(names(lakes73clean) %in% c("RR","RRn","BdgTn"))]

pearson_pooled <- predvar %>% as.data.frame() %>% setNames("param")
estimates <- c()
p_values <- c()
for(y in 1:length(predvar)){
  cor_formula <- as.formula(paste("~RR+",predvar[y],sep=""))
  for(x in 1:M){
  df <- complete(RRmice,x) %>% dplyr::select(!c("RRn","BdgTn"))
  dfcor <- cor.test(formula = cor_formula, data = df)
  estimates[x] <- dfcor$estimate[[1]]
  p_values[x] <- dfcor$p.value[[1]]
}

cor_pooled <- pool.scalar(estimates,p_values,n=73,k=1) 
pearson_pooled$r[y] <- cor_pooled$qbar
pearson_pooled$p_value[y] <- cor_pooled$ubar

}
pearson_pooled <- pearson_pooled[order(pearson_pooled$r,decreasing = T),]
g <- ggplot(data = filter(pearson_pooled,p_value <= 0.05),aes(x=reorder(param,order(r,decreasing=T)),y=r,label=param))+
  geom_col(aes(fill=p_value))+ 
  geom_text(angle=90,hjust=0,nudge_y = 0.01,size=7)+
  theme_bw(base_size=24)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.grid.minor = element_blank())+
  scale_fill_gradientn(colors = c("skyblue","dodgerblue4","firebrick"),
                       breaks = c(0,0.01,0.025, 0.05),labels = c("0","0.01","0.025","0.05"),limits = c(0,0.05), 
                       guide="legend")+
  labs(x = "", y= "r", title = "Pearson correlation coefficient with RR") + ylim(-0.5,1)
  print(g)
ggsave("8.version_control/r_RR.png",width=10,height = 6)

# correlation plot for parameters with C ----
png("C_corplot.png")
C_df <- select(lakes73,c("RR","TOC","DOC","CN","CP"))
corrplot::corrplot(corr = cor(C_df), is.corr = FALSE, p.mat = as.matrix(cor_pmat(C_df)[2:length(cor_pmat(C_df))]), method = "number",type = "upper",tl.cex = 2, tl.col = "black", number.cex = 2,hclust="ward",cl.pos = "n")
dev.off()

png("RR_s_corplot.png",width = 600, height = 600)
covar_RR <- filter(spearman_pooled,p_value <= 0.05)%>% filter(abs(r) > 0.3) %>% pull("param")
RR_scor <- micombine.cor(RRmice,variables = c("RR",covar_RR), method="spearman") 
corrplot::corrplot(corr = attr(RR_scor,"r_matrix"), p.mat = attr(RR_scor,"p_matrix"),sig.level = 0.005, order="AOE",
                   method = "number",type = "upper",tl.cex = 1.5, tl.col = "black", number.cex = 1.5,cl.pos = "n")
dev.off()

# -----
# cor.test on mice for BdgT -----

predvar <- names(lakes73clean)[-which(names(lakes73clean) %in% c("BdgT","RRn"))]

pearson_pooled_bdgT <- predvar %>% as.data.frame() %>% setNames("param")
estimates <- c()
p_values <- c()
for(y in 1:length(predvar)){
  cor_formula <- as.formula(paste("~BdgT+",predvar[y],sep=""))
  for(x in 1:M){
    df <- complete(RRmice,x)
    dfcor <- cor.test(formula = cor_formula, data = df)
    estimates[x] <- dfcor$estimate[[1]]
    p_values[x] <- dfcor$p.value[[1]]
  }
  
  cor_pooled <- pool.scalar(estimates,p_values,n=73,k=1) 
  pearson_pooled_bdgT$r[y] <- cor_pooled$qbar
  pearson_pooled_bdgT$p_value[y] <- cor_pooled$ubar
  
}

pearson_pooled_bdgT <- pearson_pooled_bdgT[order(pearson_pooled_bdgT$r,decreasing = T),]

g <- ggplot(data = filter(pearson_pooled_bdgT,p_value <= 0.05),aes(x=reorder(param,order(r,decreasing=T)),y=r,label=param))+
  geom_col(aes(fill=p_value))+ 
  geom_text(angle=90,hjust=0,nudge_y = 0.01,size=7)+
  theme_bw(base_size=28)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.grid.minor = element_blank())+
  scale_fill_gradientn(colors = c("skyblue","dodgerblue4","firebrick"),
                       breaks = c(0,0.01,0.025, 0.05),labels = c("0","0.01","0.025","0.05"),limits = c(0,0.05), 
                       guide="legend")+
  xlab("")+ ylim(-0.5,1)
print(g)
ggsave("8.version_control/pearson_pooled_bdgt.png",width = 8, height = 5)


# -----
# cor.test on mice with Spearman -----
predvar <- names(lakes73clean)[-which(names(lakes73clean) %in% c("RR","RRn"))]

spearman_pooled <- predvar %>% as.data.frame() %>% setNames("param")
estimates <- c()
p_values <- c()
for(y in 1:length(predvar)){
  cor_formula <- as.formula(paste("~RR+",predvar[y],sep=""))
  for(x in 1:M){
    df <- complete(RRmice,x)
    dfcor <- cor.test(formula = cor_formula, data = df, method = "spearman",exact=F)
    estimates[x] <- dfcor$estimate[[1]]
    p_values[x] <- dfcor$p.value[[1]]
  }
  
  cor_pooled <- pool.scalar(estimates,p_values,n=73,k=1) 
  spearman_pooled$r[y] <- cor_pooled$qbar
  spearman_pooled$p_value[y] <- cor_pooled$ubar
  
}

spearman_pooled <- spearman_pooled[order(spearman_pooled$r,decreasing = T),]

g <- ggplot(data = filter(spearman_pooled,p_value <= 0.05),aes(x=reorder(param,order(r,decreasing=T)),y=r,label=reorder(param,order(r,decreasing=T))))+
  geom_col(aes(fill=p_value))+ 
  geom_text(angle=90,hjust=0,nudge_y = 0.01,size=7)+
  theme_bw(base_size=28)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.grid.minor = element_blank())+
  scale_fill_gradientn(colors = c("skyblue","dodgerblue4","firebrick"),
                       breaks = c(0,0.01,0.025, 0.05),labels = c("0","0.01","0.025","0.05"),limits = c(0,0.05), 
                       guide="legend")+
  xlab("")+ ylim(-1,1)
print(g)
ggsave("8.version_control/spearman_RR.png",width=8,height = 5)

# -----
# cor.test Spearman on mice for BdgT -----

predvar <- names(lakes73clean)[-which(names(lakes73clean) %in% c("BdgT","RR"))]

spearman_pooled_bdgT <- predvar %>% as.data.frame() %>% setNames("param")
estimates <- c()
p_values <- c()
for(y in 1:length(predvar)){
  cor_formula <- as.formula(paste("~BdgT+",predvar[y],sep=""))
  for(x in 1:M){
    df <- complete(RRmice, x)
    dfcor <- cor.test(formula = cor_formula, data = df,method = "spearman",exact = F)
    estimates[x] <- dfcor$estimate[[1]]
    p_values[x] <- dfcor$p.value[[1]]
  }
  
  cor_pooled <- pool.scalar(estimates,p_values,n=73,k=1) 
  spearman_pooled_bdgT$r[y] <- cor_pooled$qbar
  spearman_pooled_bdgT$p_value[y] <- cor_pooled$ubar
  
}

spearman_pooled_bdgT <- spearman_pooled_bdgT[order(spearman_pooled_bdgT$r,decreasing = T),]

g <- ggplot(data = filter(spearman_pooled_bdgT,p_value <= 0.05),aes(x=reorder(param,order(r,decreasing=T)),y=r,label=reorder(param,order(r,decreasing=T))))+
  geom_col(aes(fill=p_value))+ 
  geom_text(angle=90,hjust=0,nudge_y = 0.01,size=7)+
  theme_bw(base_size=28)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.grid.minor = element_blank())+
  scale_fill_gradientn(colors = c("skyblue","dodgerblue4","firebrick"),
                       breaks = c(0,0.01,0.025, 0.05),labels = c("0","0.01","0.025","0.05"),limits = c(0,0.05), 
                       guide="legend")+
  xlab("")+ ylim(-1,1)
print(g)
ggsave("8.version_control/spearman_bdgt.png",width=8,height = 5)

# mice cor BdgT ----
png("BdgT_cor.png")
corvar_BdgT <- filter(pearson_pooled_bdgT,p_value <=0.05) %>% pull("param")
BdgT_cor <- micombine.cor(RRmice,variables = c("BdgT",corvar_BdgT))
corrplot::corrplot(corr = attr(BdgT_cor,"r_matrix"),  p.mat = attr(BdgT_cor,"p_value"), method = "number",type = "upper",tl.cex = 1.5, tl.col = "black", number.cex = 1.5,hclust="ward",cl.pos = "n")
dev.off()

png("BdgT_s_corplot.png",width = 600, height = 600)
covars_BdgT <- filter(spearman_pooled_bdgT,p_value <= 0.05) %>% filter(abs(r) > 0.3)  %>% pull("param")
BdgT_scor <- micombine.cor(RRmice,variables = c("BdgT",covars_BdgT), method="spearman") 
corrplot::corrplot(corr = attr(BdgT_scor,"r_matrix"), p.mat = attr(BdgT_scor,"p_matrix"),sig.level = 0.005, order="AOE",
                   method = "number",type = "upper",tl.cex = 1.5, tl.col = "black", number.cex = 1.5,cl.pos = "n")
dev.off()
# -----
# lasso function for mice objects -----
milasso <- function(mice.object,M,resp.var,excluded.var = "none"){
  
  n <- dim(mice.object$data)[1] 
  lasso_coef <- c("Intercept",names(select(mice.object$data,!c(resp.var,excluded.var)))) %>% as.data.frame() %>% setNames("param")
  lasso_pred <- c(rep(NA,n)) %>% as.data.frame()
  
  
  for (z in 1:M){
    lasso_coef %>% tibble::add_column(z=NA)
    data <- complete(mice.object,z) %>% select(!excluded.var)
    set.seed(5)
    lasso.formula <- as.formula(paste(resp.var,"~.",sep=""))
    fit <- cv.glmnet(formula = lasso.formula,data = data)
    coef_fit <- fit %>% coef(s="lambda.min") %>% summary()
    
    for (w in 1:dim(coef_fit)[1]){
      lasso_coef[coef_fit[w,1],z+1] <- coef_fit[w,3]
    }
    
    newmat <- complete(mice.object,z) %>% select(!c(resp.var,excluded.var)) 
    lasso_pred[,z] <- predict(fit,newmat[1:n,],s="lambda.min")
  }
  
  lasso_coef$pooled <- NA
  lasso_coef$number <- NA
  lasso_pred$pooled <- rowMeans(lasso_pred)
  
  
  for(x in 1:length(lasso_coef$param)){
    estimates <- lasso_coef[x,] %>% dplyr::select(c(1:M+1)) %>% as.numeric()
    lasso_coef$pooled[x] <- pool.scalar(estimates,c(rep(1,M)),n=n,k=1)$qbar 
    lasso_coef$number[x] <- M - sum(is.na(estimates))
  }

  final <- list(lasso_coef,lasso_pred)
  names(final) <- c("lasso_coef","lasso_pred")
  return(final)
  
}

# lasso calculations -----
lasso_RR <- milasso(RRmice,M,"RR","RRn")
lasso_BdgT <- milasso(RRmice,M,"BdgT","RRn")
lasso_RRn <- milasso(RRmice,M,"RRn",c("RR","DOC"))

rmse(lakes73$RR,lasso_RR$lasso_pred$pooled)
rmse(lakes73$BdgT,lasso_BdgT$lasso_pred$pooled)
rmse(lakes73clean$RRn,lasso_RRn$lasso_pred$pooled)

write_xlsx(lasso_RR$lasso_coef,"8.version_control/lasso_coef_RR.xlsx")
write_xlsx(lasso_BdgT$lasso_coef,"8.version_control/lasso_coef_BdgT.xlsx")
write_xlsx(lasso_RRn$lasso_coef,"8.version_control/lasso_coef_RRn.xlsx")

# -----
# lasso on log -----
lakes73log <- log(lakes73clean)
lakes73log[lakes73log == "-Inf"] <- NA
M <- 10
set.seed(5)
RRlogmice <- lakes73log %>% mice(method = "cart", m = M)


lasso_RRlog <- milasso(RRlogmice,M,"RR","BdgT")
write_xlsx(lasso_RRlog$lasso_coef,"8.version_control/lasso_RRlog.xlsx")

plot(lakes73log$RR,lasso_RRlog$lasso_pred$pooled)+abline(a=0,b=1)
rmse(lakes73log$RR,lasso_RRlog$lasso_pred$pooled)
lm(data=lakes73log,RR~Na*Mg*Fe+SUVA*SARuv*SARvis+CH4+CN*CP) %>% summary()

lasso_BdgTlog <- milasso(RRlogmice,M,"BdgT","RR")
write_xlsx(lasso_BdgTlog$lasso_coef, "8.version_control/lasso_BdgTlog.xlsx")

plot(lakes73log$BdgT,lasso_BdgTlog$lasso_pred$pooled)+abline(a=0,b=1)
rmse(lakes73log$BdgT,lasso_BdgTlog$lasso_pred$pooled)
lm(data=lakes73log,BdgT~Cells+TOC+DN+pH+Alkalinity+Ca+Fe+Cd+Mn+SARuv+SARvis+O2) %>% summary()


# ----- 

# Gauss-lasso estimator RR -----
lm(RR~CN,data=lakes73) %>% summary()
ggplot(lakes73)+geom_text(aes(x=RR,y=1.09+0.083*CN,label=Lake_ID))+geom_abline(slope=1,intercept=0)

# -----
# Gauss-lasso estimator BdgT -----
lm(BdgT~DN+pH,data=lakes73clean) %>% summary()
ggplot(lakes73)+geom_text(aes(x=BdgT,y=-3.46+5.22*DN+0.83*pH,label=Lake_ID))+geom_abline(slope=1,intercept=0)

lm(BdgT~pH+DN+Cells,data=lakes73) %>% summary()
lm(BdgT~DN,lakes73) %>% summary()

cor.test(formula = BdgT ~ 2.73+2.82*DN+0.15*pH, data = lakes73clean)

# cor.test on mice for RR/DOC -----

lakes73clean$RRn <- lakes73clean$RR/lakes73clean$DOC
predvar <- names(lakes73clean)[-which(names(lakes73clean) %in% c("BdgT","RR","RRn"))]


pearson_pooled_RRn <- predvar %>% as.data.frame() %>% setNames("param")
estimates <- c()
p_values <- c()
for(y in 1:length(predvar)){
  cor_formula <- as.formula(paste("~RRn+",predvar[y],sep=""))
  for(x in 1:M){
    df <- complete(RRmice,x)
    dfcor <- cor.test(formula = cor_formula, data = df)
    estimates[x] <- dfcor$estimate[[1]]
    p_values[x] <- dfcor$p.value[[1]]
  }
  
  cor_pooled <- pool.scalar(estimates,p_values,n=73,k=1) 
  pearson_pooled_RRn$r[y] <- cor_pooled$qbar
  pearson_pooled_RRn$p_value[y] <- cor_pooled$ubar
  
}
pearson_pooled_RRn <- pearson_pooled_RRn[order(pearson_pooled_RRn$r,decreasing = T),]

g <- ggplot(data = filter(pearson_pooled_RRn,p_value <= 0.04),aes(x=reorder(param,order(r,decreasing=T)),y=r,label=param))+
  geom_col(aes(fill=p_value))+ 
  geom_text(angle=90,hjust=0,nudge_y = 0.01,size=7)+
  theme_bw(base_size=28)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.grid.minor = element_blank())+
  scale_fill_gradientn(colors = c("skyblue","dodgerblue4","firebrick"),
                       breaks = c(0,0.02,0.04, 0.06),labels = c("0","0.02","0.04","0.06"),limits = c(0,0.06), 
                       guide="legend")+
  xlab("")+ ylim(-0.5,0.5)
print(g)
ggsave("8.version_control/pearson_RRn.png", width=8,height = 5)

spearman_pooled_RRn <- predvar %>% as.data.frame() %>% setNames("param")
estimates <- c()
p_values <- c()
for(y in 1:length(predvar)){
  cor_formula <- as.formula(paste("~RRn+",predvar[y],sep=""))
  for(x in 1:M){
    df <- complete(RRmice,x)
    dfcor <- cor.test(formula = cor_formula, data = df, method = "spearman", exact = FALSE)
    estimates[x] <- dfcor$estimate[[1]]
    p_values[x] <- dfcor$p.value[[1]]
  }
  
  cor_pooled <- pool.scalar(estimates,p_values,n=73,k=1) 
  spearman_pooled_RRn$r[y] <- cor_pooled$qbar
  spearman_pooled_RRn$p_value[y] <- cor_pooled$ubar
  
}

spearman_pooled_RRn <- spearman_pooled_RRn[order(spearman_pooled_RRn$r,decreasing = T),]
g <- ggplot(data = filter(spearman_pooled_RRn,p_value <= 0.05),aes(x=reorder(param,order(r,decreasing=T)),y=r,label=param))+
  geom_col(aes(fill=p_value))+ 
  geom_text(angle=90,hjust=0,nudge_y = 0.01,size=7)+
  theme_bw(base_size=28)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.grid.minor = element_blank())+
  scale_fill_gradientn(colors = c("skyblue","dodgerblue4","firebrick"),
                       breaks = c(0,0.01,0.025, 0.05),labels = c("0","0.01","0.025","0.05"),limits = c(0,0.05), 
                       guide="legend")+
  xlab("")+ ylim(-1,1)
print(g)
ggsave("8.version_control/spearman_RRn.png",width=15,height = 5)

# -----
# miceadds corplots RRn -----
png("RRn_corplot.png")
RRn_cor <- micombine.cor(RRmice,variables = c("RRn","DOC","SUVA","SR","sVISa","DN","As","V","RR")) 
corrplot::corrplot(corr = attr(RRn_cor,"r_matrix"), p.mat = attr(RRn_cor,"p_matrix"), order="hclust",hclust="centroid",
                   method = "number",type = "upper",tl.cex = 1.5, tl.col = "black", number.cex = 1.5,cl.pos = "n")
dev.off()

png("RRn_s_corplot.png",width = 600, height = 600)
covar_RRn <- filter(spearman_pooled_RRn,p_value <= 0.05) %>% filter(abs(r) > 0.5) %>% pull("param")
RRn_scor <- micombine.cor(RRmice,variables = c("RRn","DOC","RR",covar_RRn), method="spearman") 
corrplot::corrplot(corr = attr(RRn_scor,"r_matrix"), p.mat = attr(RRn_scor,"p_matrix"),sig.level = 0.005, order="hclust",hclust="ward.D2",
                   method = "number",type = "upper",tl.cex = 1.5, tl.col = "black", number.cex = 1.5,cl.pos = "n")
dev.off()

# -----

# Gauss - lasso estimates RRn-----
lm(RRn~DN+SUVA+SR+V,data=lakes73clean) %>% summary()
lm(RRn~DN+SR,data=lakes73clean) %>% summary()

# -----
