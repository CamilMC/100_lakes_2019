#-----
# Load libraries & functions
#-----


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
library(Metrics)
library(rlist)
library(plotmo)

source("8.version_control/100lakes_analysis_functions.R")
source("8.version_control/bacteria.R")

norgemap <- st_read("3.Maps/TM_WORLD_BORDERS-0.3.shp")

library("rnaturalearth")
library(rgeos)
world <- ne_countries(scale = "medium", returnclass = "sf")

#-----
# Load data
#-----

lakes73 <- read.csv("8.version_control/rr100lakes.csv")
lakes73$SARvis[which(lakes73$SARvis == "Inf")] <- "NA"
lakes73num <- lakes73 %>% sapply(as.numeric) %>% as.data.frame() %>% select_if(not_all_na)

#-----
# Descriptive stats
#-----
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
       
i <- ggplot()+geom_boxplot(data=lakes73,aes(x="",y=log(RR)))
outliers <- ggplot_build(i)[["data"]][[1]][["outliers"]] %>% as.data.frame()
outliers_ID <-  filter(lakes73, log(lakes73$RR) %in% outliers[[1]]) %>% dplyr::select(Lake_name)
outliers$Lake_ID <- outliers_ID$Lake_name
names(outliers) <- c("RR","Lake_name")
i <- ggplot()+geom_boxplot(data=lakes73,aes(x="",y=log(RR)),outlier.size = 3, width = 1.5)+theme_light(base_size = 20)+xlab("")+ylab("log(RR)")+
  geom_text_repel(data=outliers,aes(x="",y=outliers[,"RR"],label=Lake_name),position = position_dodge(width=2),size = 6, col = "steelblue4",segment.size = 1)
print(i)
ggsave("boxplot_logRR.png",width=4,height=6)

j <- ggplot()+geom_boxplot(data=lakes73,aes(x="",y=log(BdgT)))
outliers <- ggplot_build(j)[["data"]][[1]][["outliers"]] %>% as.data.frame()
outliers_ID <-  filter(lakes73, log(lakes73$BdgT) %in% outliers[[1]]) %>% dplyr::select(Lake_name)
outliers$Lake_ID <- outliers_ID$Lake_name
names(outliers) <- c("BdgT","Lake_name")
j <- ggplot()+geom_boxplot(data=lakes73,aes(x="",y=log(BdgT)),outlier.size = 3, width = 1.5)+theme_light(base_size = 20)+xlab("")+ylab("log(BdgT)")+
  geom_text_repel(data=outliers,aes(x="",y=outliers[,"BdgT"],label=Lake_name),position = position_dodge(width=2),size = 6, col = "steelblue4",segment.size = 1)
print(j)
ggsave("boxplot_logBdgT.png",width=4,height=6)


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


# distribution of data ----
pdf("8.version_control/distribution.pdf")
hist.df(lakes73num)
dev.off()

# corrplot -----

corlakes73 <- cor(lakes73num, use="complete.obs",method = c("pearson"))
write_xlsx(as.data.frame(corlakes73),"8.version_control/corlakes73.xlsx")
corlakes73_pvalues <- corrplot::cor.mtest(lakes73num,conf.level=0.95)
corrplot::corrplot(corlakes73,method = c("number"),type=c("upper"),tl.cex = 0.6,number.cex=0.45,p.mat = corlakes73_pvalues$p,insig = "blank")

print.cor.signif2(lakes73,"RR")

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

# correLations ----
pdf("5.100_lakes/all.parameters.correlations.pdf",width=15,height = 8)
sapply(names(lakes73),print.cor.signif2, df=lakes73)
dev.off()
# -----
# Maps
# -----

ggplot(lakes73)+geom_text(aes(x=incub_date,y=RR,col=incub_date,label=Lake_ID),size=5,show.legend = F)+
  theme_light(base_size=20)+theme(legend.position = "none",axis.text.x = element_text(angle = 45))+
  scale_color_brewer(palette="Paired")
ggplot(lakes73)+geom_text(aes(x=incub_date,y=RR,col=DOC,label=Lake_ID),size=5,show.legend = T)+
  theme_light(base_size=20)+theme(axis.text.x = element_text(angle = 30))+
  scale_color_viridis_c(direction = -1, end = 0.8)
ggplot() + geom_sf(data=world,fill = "white") + coord_sf(xlim = c(4.5,12.9), ylim = c(58, 62), expand = F)+
  geom_point(data=lakes73,aes(x=Long,y=Lat,col=incub_date,size=RR))+
  geom_point(data=lakes8,aes(x=Long,y=Lat),col="black",size=5,shape="*")+
  scale_color_brewer(palette="Paired")+
  theme_void(base_size = 24)+labs(col="Incubation date")
ggplot() + geom_sf(data=world,fill = "white") + coord_sf(xlim = c(4.5,12.9), ylim = c(58, 62), expand = F)+
  geom_point(data=lakes73,aes(x=Long,y=Lat,col=incub_date,size=DOC))+
  geom_point(data=lakes8,aes(x=Long,y=Lat),col="black",size=5,shape="*")+
  scale_color_brewer(palette="Paired")+
  theme_void(base_size = 24)+labs(col="Incubation date")
ggplot(lakes73num)+geom_text(aes(x=BdgT,y=RR,label=Lake_ID))+theme_light(base_size=20)
ggplot(lakes73num)+geom_text(aes(x=SARuv,y=BdgT,label=Lake_ID))+theme_light(base_size=20)


pdf("5.100_lakes/RR_maps.pdf",width=10,height=10)

ggplot()+geom_sf(data=norgemap)+ geom_point(data=lakes73,aes(x=Long,y=Lat,col=RR),size=5)+ 
  geom_point(data=filter(lakes73,RR >8),aes(x=Long,y=Lat,size=RR))+scale_size(limits=c(8,30),breaks=c(10,20,30),range=c(4,6))+
  theme_void(base_size = 24)+labs(x="",y="",col=expression(atop("RR",paste("(",mu,"mol/L/h)"))),size="outliers")+
  scale_color_gradientn(colors = viridis(9,direction = -1,end=0.9),limits=c(0,8))+
  #geom_text(data=filter(lakes73,RR>8),aes(x=Long,y=Lat,label=round(RR,0)),nudge_y=0.1)+
  guides(color=guide_legend(order=1),size = guide_legend(order=2))+
  theme(legend.position = c(.87,.45))+
  xlim(5,15)+ylim(58,63)
ggplot()+geom_sf(data=norgemap)+ geom_point(data=lakes73,aes(x=Long,y=Lat,col=RR.OD),size=5)+
  theme_void(base_size = 24)+labs(x="",y="",col=expression(atop("RR/OD",(h^{-1}))))+
  scale_color_viridis_c(direction = -1, end = 0.9)+
  theme(legend.position = c(.9,.45))+
  xlim(5,15)+ylim(58,63)

ggplot()+geom_sf(data=norgemap)+ geom_point(data=lakes73,aes(x=Long,y=Lat,col=RR.OD),size=5)+
  geom_point(data=filter(lakes73,RR.OD >2),aes(x=Long,y=Lat,size=RR.OD))+scale_size(limits=c(2,10),breaks=c(3,5,10),range=c(4,6))+
  theme_void(base_size = 24)+labs(x="",y="",col=expression(atop("RR/DOC",(h^{-1}))),size="outliers")+
  scale_color_gradientn(colors = viridis(6,direction = -1,end=0.9),limits=c(0,2))+
  #geom_text(data=filter(lakes73,RR.OD >2),aes(x=Long,y=Lat,label=round(RR.OD,2)),nudge_y = -0.09,nudge_x = -0.1)+
  guides(color=guide_legend(order=1),size = guide_legend(order=2))+
  theme(legend.position = c(.87,.45))+
  xlim(5,15)+ylim(58,63)

ggplot()+geom_sf(data=norgemap)+ geom_point(data=lakes73,aes(x=Long,y=Lat,col=OD.DOC),size=4)+
  theme_void(base_size = 24)+labs(x="",y="",col="%DOC\nconsumed",size="outliers")+
  scale_color_gradientn(colors = viridis(6,direction = -1,end=0.9),limits=c(0,4))+
  geom_point(data=filter(lakes73,OD.DOC >4),aes(x=Long,y=Lat,size=OD.DOC))+scale_size(limits=c(4,20),breaks=c(5,10,15),range=c(4,6))+
  #  geom_text(data=filter(lakes73,OD.DOC >3),aes(x=Long,y=Lat,label=round(OD.DOC,2)),nudge_y = -0.09,nudge_x = -0.1)+  
  guides(color=guide_legend(order=1),size = guide_legend(order=2))+
  theme(legend.position = c(.87,.5))+
  xlim(5,15)+ylim(58,63)

ggplot()+geom_sf(data=norgemap)+ geom_point(data=lakes73,aes(x=Long,y=Lat,col=OD),size=4)+
  theme_void(base_size = 24)+labs(x="",y="",col="OD\n(umol/L)")+
  scale_color_gradientn(colors = viridis(6,direction = -1,end=0.9),limits=c(0,30))+
  geom_text(data=filter(lakes73,OD >30),aes(x=Long,y=Lat,label=round(OD,0)),nudge_y = -0.09,nudge_x = -0.1)+  
  theme(legend.position = c(.9,.7))+
  xlim(5,15)+ylim(58,63)

ggplot()+geom_sf(data=norgemap)+ geom_point(data=lakes73,aes(x=Long,y=Lat,col=DN),size=4)+
  geom_point(data=filter(lakes73,DN > 0.6),aes(x=Long,y=Lat,size= DN))+scale_size(limits=c(0.6,0.9))+
  theme_void(base_size = 24)+labs(x="",y="",col="DN\nmg/L",size="Outliers")+
  scale_color_gradientn(colors = viridis(6,direction = -1,end=0.9),limits=c(0,0.6))+
  #  geom_text(data=filter(lakes73,DN > 0.6),aes(x=Long,y=Lat,label=round(DN,1)),nudge_y = -0.09,nudge_x = -0.1)+  
  theme(legend.position = c(.9,.6))+
  xlim(5,15)+ylim(58,63)

ggplot()+geom_sf(data=norgemap)+ geom_point(data=lakes73,aes(x=Long,y=Lat,col=DP),size=4)+
  geom_point(data=filter(lakes73,DP > 10),aes(x=Long,y=Lat,size= DP))+scale_size(limits=c(10,16))+
  theme_void(base_size = 24)+labs(x="",y="",col="DP\nug/L",size="Outliers")+
  scale_color_gradientn(colors = viridis(6,direction = -1,end=0.9),limits=c(0,10))+
  #  geom_text(data=filter(lakes73,DN > 0.6),aes(x=Long,y=Lat,label=round(DN,1)),nudge_y = -0.09,nudge_x = -0.1)+  
  theme(legend.position = c(.9,.6))+
  xlim(5,15)+ylim(58,63)

ggplot() + geom_sf(data=world,fill = "white") + coord_sf(xlim = c(4.5,12.9), ylim = c(58, 62), expand = F)+
  geom_point(data=lakes73,aes(x=Long,y=Lat,col=CN),size=3)+
  geom_point(data=filter(lakes73,CN > 0.050),aes(x=Long,y=Lat,size= CN))+scale_size(limits=c(0.05,0.5))+
  scale_color_gradientn(colors = viridis(6,direction = -1,end=0.9),limits=c(0,0.05))+
  theme_void(base_size = 20)+labs(col="CN")

ggplot()+geom_sf(data=norgemap)+ geom_point(data=lakes73,aes(x=Long,y=Lat,col=CN),size=4)+
  #  geom_point(data=filter(lakes73,CN > 50),aes(x=Long,y=Lat,size= CN))+scale_size(limits=c(50,150),range=c(4,6),breaks=c(50,100,150))+
  theme_void(base_size = 24)+labs(x="",y="",col="CN",size="Outliers")+
  # scale_color_gradientn(colors = viridis(6,direction = -1,end=0.9))+
  scale_color_gradientn(colors = viridis(6,direction = -1,end=0.9),limits=c(8,50))+
  # geom_text(data=filter(lakes73,DN > 0.6),aes(x=Long,y=Lat,label=round(DN,1)),nudge_y = -0.09,nudge_x = -0.1)+ 
  guides(color=guide_legend(order=1),size = guide_legend(order=2))+
  theme(legend.position = c(.9,.6))+
  xlim(5,15)+ylim(58,63)

ggplot()+geom_sf(data=norgemap)+ geom_point(data=lakes73,aes(x=Long,y=Lat,col=DOC),size=5)+
  #geom_point(data=filter(lakes73,DOC > 20),aes(x=Long,y=Lat,size= DOC))+scale_size(limits=c(20,21))+
  theme_void(base_size = 26)+labs(x="",y="",col="DOC\nmg/L",size="Outliers")+
  # scale_color_gradientn(colors = viridis(6,direction = -1,end=0.9))+
  scale_color_gradientn(colors = viridis(6,direction = -1,end=0.9),limits=c(1,20))+
  #  geom_text(data=filter(lakes73,DOC > 15),aes(x=Long,y=Lat,label=round(DOC,1)),nudge_y = -0.1,nudge_x = -0.0)+ 
  guides(color=guide_legend(order=1),size = guide_legend(order=2))+
  theme(legend.position = c(.75,.2))+
  xlim(5,15)+ylim(58,63)

dev.off()


ggplot() + geom_sf(data=world,fill = "white") + coord_sf(xlim = c(4.5,12.9), ylim = c(58, 62), expand = F)+
  geom_point(data=lakes73,aes(x=Long,y=Lat,col=Lake_ID),size=3)+
  geom_text_repel(data=lakes73,aes(x=Long,y=Lat,col=Lake_ID,label=Lake_ID),nudge_y=0.1)+
  scale_color_viridis_c()+theme_void(base_size = 20)

ggsave("8.version_control/map_LakeID.jpeg",width=15,height=10)

ggplot() + geom_sf(data=world,fill = "white") + coord_sf(xlim = c(4.5,12.9), ylim = c(58, 62), expand = F)+
  geom_point(data=lakes73,aes(x=Long,y=Lat,col=DOC,size=DOC))+
  scale_color_gradientn(colors = viridis(6,direction = -1, begin = 0.3, end=0.9))+scale_size(range = c(2.5,8))+
  guides(color = guide_legend(),size=guide_legend())+theme_void(base_size = 24)
ggsave("8.version_control/map_DOC.jpeg",width=15,height=10)

#-----
# MICE on lakes73 
#-----
# create multiple imputation datasets to fill NA -----
param2remove <- c("auc","width","OD","Lake_ID","NVE_number","NIVA_date","X","Long","Lat","Altitude","CBA_week","tmax","tmax_h","ox_initial","F","pH_bio","EC_bio","NO2","NO3","DNA","p_N2","p_O2","c_O2","p_CO2","c_CO2","p_CH4","p_N2O","d18O","d2H","Tso","a254","a400","a410","a600","a275","a295","s_275_295","a350","a400","s_350_400","lag_bdg","H","NP","CNP","dist_bact","PO4")
param2keep <- c("RR","BdgT","T","pH","EC","Cells","DOC","DN","DP","CN","CP","Alkalinity","Cl","B","SO4","Na","K","Ca","Mg","Fe","Al","SUVA","SARuv","SARvis","SR","sVISa","O2","CO2","CH4","N2O","p_O2","c_O2","p_CO2","c_CO2","p_CH4")
lakes73u <- lakes73num %>% dplyr::select(which(names(lakes73num)%in% param2keep))
lakes73u$RRn <- lakes73u$RR/lakes73u$DOC
write_xlsx(lakes73u,"8.version_control/lakes73u.xlsx")

M <- 50
set.seed(5)
RRmice <- lakes73u %>% mice(method = "cart",m = M)
# lasso function for mice objects -----
milasso <- function(mice.object,M,resp.var,excluded.var){
  
  n <- dim(mice.object$data)[1] 
  lasso_coef <- c("Intercept",names(dplyr::select(mice.object$data,!c(resp.var,excluded.var)))) %>% as.data.frame() %>% setNames("param")
  lasso_pred <- c(rep(NA,n)) %>% as.data.frame()
  
  
  for (z in 1:M){
    lasso_coef %>% tibble::add_column(z=NA)
    data <- complete(mice.object,z) %>% dplyr::select(!excluded.var)
    set.seed(5)
    lasso.formula <- as.formula(paste(resp.var,"~.",sep=""))
    fit <- cv.glmnet(formula = lasso.formula,data = data)
    coef_fit <- fit %>% coef(s="lambda.min") %>% summary()
    
    for (w in 1:dim(coef_fit)[1]){
      lasso_coef[coef_fit[w,1],z+1] <- coef_fit[w,3]
    }
    
    newmat <- complete(mice.object,z) %>% dplyr::select(!c(resp.var,excluded.var)) 
    lasso_pred[,z] <- predict(fit,newmat[1:n,],s="lambda.min")
  }
  
  lasso_coef$pooled <- NA
  lasso_coef$number <- NA
  lasso_pred$pooled <- rowMeans(lasso_pred)
  
  
  for(x in 1:length(lasso_coef$param)){
    estimates <- lasso_coef[x,] %>% dplyr::select(c(1:M+1)) %>% as.numeric()
    lasso_coef$pooled[x] <- mean(estimates,na.rm = T)
    lasso_coef$number[x] <- M - sum(is.na(estimates))
  }

  final <- list(lasso_coef,lasso_pred)
  names(final) <- c("lasso_coef","lasso_pred")
  return(final)
  
}

milasso2 <- function(mice.object,M,resp.var,excluded.var){
  
  n <- dim(mice.object$data)[1] 
  lasso_coef <- c("Intercept",names(dplyr::select(mice.object$data,!c(resp.var,excluded.var)))) %>% as.data.frame() %>% setNames("param")
  lasso_pred <- c(rep(NA,n)) %>% as.data.frame()
  
  
  for (z in 1:M){
    lasso_coef %>% tibble::add_column(z=NA)
    data <- complete(mice.object,z) %>% dplyr::select(!excluded.var)
    set.seed(5)
    lasso.formula <- as.formula(paste(resp.var,"~.",sep=""))
    fit <- cv.glmnet(formula = lasso.formula,data = data)
    coef_fit <- fit %>% coef(s="lambda.min") %>% summary()
    
    for (w in 1:dim(coef_fit)[1]){
      lasso_coef[coef_fit[w,1],z+1] <- coef_fit[w,3]
    }
    
    newmat <- complete(mice.object,z) %>% dplyr::select(!c(resp.var,excluded.var)) 
    lasso_pred[,z] <- predict(fit,newmat[1:n,],s="lambda.min")
  }
  
  lasso_coef$pooled <- NA
  lasso_coef$number <- NA
  lasso_pred$pooled <- rowMeans(lasso_pred)
  
  
  for(x in 1:length(lasso_coef$param)){
    estimates <- lasso_coef[x,] %>% dplyr::select(c(1:M+1)) %>% as.numeric()
    lasso_coef$pooled[x] <- mean(estimates,na.rm = T)
    lasso_coef$number[x] <- M - sum(is.na(estimates))
  }
  
  final <- list(lasso_coef,lasso_pred)
  names(final) <- c("lasso_coef","lasso_pred")
  return(final)
  
}

#-----
# Pearson for RR -----
predvar <- names(lakes73u)[-which(names(lakes73u) %in% c("RR","RRn"))]

pearson_pooled <- predvar %>% as.data.frame() %>% setNames("param")
estimates <- c()
p_values <- c()
for(y in 1:length(predvar)){
  cor_formula <- as.formula(paste("~RR+",predvar[y],sep=""))
  for(x in 1:M){
    df <- complete(RRmice,x) %>% dplyr::select(!c("RRn"))
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

png("8.version_control/RR_corplot.png")
RR_corvar <- filter(pearson_pooled,p_value <= 0.05) %>% pull("param")
RR_cor <- micombine.cor(RRmice,variables = c("RR",RR_corvar), method="pearson") 
matrix_RR_cor <- attr(RR_cor,"r_matrix")
colnames(matrix_RR_cor) <- c("DOC","RR", "C:N","C:P")
rownames(matrix_RR_cor) <- c("DOC","RR", "C:N","C:P")
matrix_RR_p <-attr(RR_cor,"p_value")
corrplot::corrplot(corr = matrix_RR_cor,  p.mat = matrix_RR_p, method = "number",type = "upper",tl.cex = 2, tl.col = "black", number.cex = 2,hclust="ward",cl.pos = "n")
dev.off()


# Spearman for RR ----

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

png("8.version_control/RR_s_corplot.png",width = 600, height = 600)
covar_RR <- filter(spearman_pooled,p_value <= 0.05) %>% pull("param")
RR_scor <- micombine.cor(RRmice,variables = c("RR",covar_RR), method="spearman") 
corrplot::corrplot(corr = attr(RR_scor,"r_matrix"), p.mat = attr(RR_scor,"p_matrix"),sig.level = 0.005, order="AOE",
                   method = "number",type = "upper",tl.cex = 1.5, tl.col = "black", number.cex = 1.5,cl.pos = "n")
dev.off()




# lasso RR  -----
lasso_RR <- milasso(RRmice,M,"RR","RRn")
rmse(lakes73$RR,lasso_RR$lasso_pred$pooled)
write_xlsx(lasso_RR$lasso_coef,"8.version_control/lasso_coef_RR.xlsx")

# Gauss-lasso estimates RR -----
lm_RR <- with(RRmice,lm(RR~CN))
lm_RR_pooled <- pool(lm_RR) %>% summary()
lm_RR_predict <- with(RRmice,predict(lm(RR~CN)))
RR_predicted <- list.cbind(lm_RR_predict$analyses) %>% rowMeans()
rmse(lakes73$RR,RR_predicted)


# Pearson for BdgT -----

predvar <- names(lakes73u)[-which(names(lakes73u) %in% c("BdgT","RRn"))]

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


# Spearman for BdgT -----

predvar <- names(lakes73u)[-which(names(lakes73u) %in% c("BdgT","RR"))]

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
png("8.version_control/BdgT_cor.png")
corvar_BdgT <- filter(pearson_pooled_bdgT,p_value <=0.05) %>% pull("param")
BdgT_cor <- micombine.cor(RRmice,variables = c("BdgT",corvar_BdgT))
corrplot::corrplot(corr = attr(BdgT_cor,"r_matrix"),  p.mat = attr(BdgT_cor,"p_value"), method = "number",type = "upper",tl.cex = 1.3, tl.col = "black", number.cex = 1.3,hclust="ward",cl.pos = "n")
dev.off()

png("8.version_control/BdgT_s_corplot.png",width = 600, height = 600)
covars_BdgT <- filter(spearman_pooled_bdgT,p_value <= 0.05) %>% pull("param")
BdgT_scor <- micombine.cor(RRmice,variables = c("BdgT",covars_BdgT), method="spearman") 
corrplot::corrplot(corr = attr(BdgT_scor,"r_matrix"), p.mat = attr(BdgT_scor,"p_matrix"),sig.level = 0.005, order="AOE",
                   method = "number",type = "upper",tl.cex = 1.5, tl.col = "black", number.cex = 1.5,cl.pos = "n")
dev.off()


# lasso BdgT -----

lasso_BdgT <- milasso(RRmice,M,"BdgT","RRn")
rmse(lakes73$BdgT,lasso_BdgT$lasso_pred$pooled)
write_xlsx(lasso_BdgT$lasso_coef,"8.version_control/lasso_coef_BdgT.xlsx")

# Gauss-lasso estimates BdgT -----
lm_BdgT <- with(RRmice,lm(BdgT~Cells+DN+pH+c_O2))
lm_BdgT_pooled <- pool(lm_BdgT) %>% summary()
lm_BdgT_predict <- with(RRmice,predict(lm(BdgT~Cells+DN+pH+c_O2)))
BdgT_predicted <- list.cbind(lm_BdgT_predict$analyses) %>% rowMeans()
rmse(lakes73$BdgT,BdgT_predicted)

# Pearson for RRn -----

predvar <- names(lakes73u)[-which(names(lakes73u) %in% c("BdgT","RRn"))]

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
  xlab("")+ ylim(-0.5,1)
print(g)
ggsave("8.version_control/pearson_RRn.png", width=8,height = 5)

# Spearman for RRn -----
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

# miceadds corplots RRn -----
png("8.version_control/RRn_corplot.png")
RRn_corvar <- filter(pearson_pooled_RRn,p_value <= 0.05) %>% pull("param")
RRn_cor <- micombine.cor(RRmice,variables = c("RRn",RRn_corvar)) 
corrplot::corrplot(corr = attr(RRn_cor,"r_matrix"), p.mat = attr(RRn_cor,"p_matrix"), order="hclust",hclust="centroid",
                   method = "number",type = "upper",tl.cex = 1.5, tl.col = "black", number.cex = 1.5,cl.pos = "n")
dev.off()

png("8.version_control/RRn_s_corplot.png",width = 600, height = 600)
covar_RRn <- filter(spearman_pooled_RRn,p_value <= 0.05) %>% pull("param")
RRn_scor <- micombine.cor(RRmice,variables = c("RRn","DOC","RR",covar_RRn), method="spearman") 
corrplot::corrplot(corr = attr(RRn_scor,"r_matrix"), p.mat = attr(RRn_scor,"p_matrix"),sig.level = 0.005, order="hclust",hclust="ward.D2",
                   method = "number",type = "upper",tl.cex = 1.5, tl.col = "black", number.cex = 1.5,cl.pos = "n")
dev.off()

# lasso RRn -----
lasso_RRn <- milasso(RRmice,M,"RRn",c("RR","DOC"))
rmse(lakes73u$RRn,lasso_RRn$lasso_pred$pooled)
write_xlsx(lasso_RRn$lasso_coef,"8.version_control/lasso_coef_RRn.xlsx")


# Gauss - lasso estimates RRn-----
lm_RRn <- with(RRmice,lm(RRn~DN+SUVA+SR))
lm_RRn_pooled <- pool(lm_RRn) %>% summary()
lm_RRn_predict <- with(RRmice,predict(lm(RRn~DN+SUVA+SR)))
RRn_predicted <- list.cbind(lm_RRn_predict$analyses) %>% rowMeans()
rmse(lakes73u$RRn,RRn_predicted)

#-----
# MICE on log lakes73 
#----- 
# mice on lakes73log -----
lakes73log <- lakes73u %>% select(!which(names(lakes73u) %in% c("pH","Cells","SUVA","SARuv","SARvis","SR","sVISa","p_O2","c_O2","p_CO2","c_CO2","p_CH4"))) %>% log10()
names(lakes73log) <- paste("log",names(lakes73log),sep="")


lakes73log <- cbind(lakes73log,select(lakes73u,c("pH","Cells","SUVA","SARuv","SARvis","SR","sVISa","p_O2","c_O2","p_CO2","c_CO2","p_CH4")))
lakes73log[lakes73log == "-Inf"] <- NA
 
pdf("8.version_control/lakes73log_distribution.pdf")
lakes73log$Lake_ID <- lakes73$Lake_ID
boxplots(lakes73log,0.95)
hist.df(lakes73log)
lakes73log$Lake_ID <- NULL
dev.off()

M <- 50
set.seed(5)

# Correlogram for all independant variables -----
indvarlogmice <- select(lakes73log,!c("logRR","logRRn","logBdgT")) %>% mice(method = "cart", m = M)

png("8.version_control/big_cor.png",width=1000,height = 1000)
lakes_cor <- micombine.cor(indvarlogmice, method="pearson") 
matrix_lakes_cor <- attr(lakes_cor,"r_matrix")
#colnames(matrix_RR_cor) <- c("log(DN)","log(DP)", "log(SO4)","log(Na)","log(K)","log(Mg)","log(RR)","log(C:N)","log(C:P)")
#rownames(matrix_RR_cor) <- c("log(DN)","log(DP)", "log(SO4)","log(Na)","log(K)","log(Mg)","log(RR)","log(C:N)","log(C:P)")
matrix_lakes_p <-attr(lakes_cor,"p_value")
corrplot::corrplot(corr = matrix_lakes_cor,  p.mat = matrix_lakes_p, method = "number",type = "upper",tl.cex = 0.8, tl.col = "black", number.cex = 0.8,hclust="ward",cl.pos = "n", sig.level = 0.99)
dev.off()

# RRlogmice -----
logcovariates <- c("logRR","logRRn","logBdgT","logDOC","logDP","logEC","logFe","logO2","logN2O","logCN","pH","Cells","SUVA","SARuv","p_O2","p_CO2","c_CO2","p_CH4")
RRlogmice <- select(lakes73log,logcovariates) %>% mice(method = "cart", m = M)

# correlogram for RRlogmice -----
png("8.version_control/correlogram.png", width= 1000, height = 1000)
lakes_cor <- micombine.cor(RRlogmice, method="pearson") 
matrix_lakes_cor <- attr(lakes_cor,"r_matrix")
colnames(matrix_lakes_cor) <- c("log(RR)","log(RRn)", "log(BdgT)","log(DOC)","log(DP)","log(DC)","log(Fe)","log(O2)","log(N2O)","log(C:N)","pH","Cells","SUVA","SARuv","p_O2","p_CO2","c_CO2","p_CH4")
rownames(matrix_lakes_cor) <- c("log(RR)","log(RRn)", "log(BdgT)","log(DOC)","log(DP)","log(DC)","log(Fe)","log(O2)","log(N2O)","log(C:N)","pH","Cells","SUVA","SARuv","p_O2","p_CO2","c_CO2","p_CH4")
matrix_lakes_p <-attr(lakes_cor,"p_value")
corrplot::corrplot(corr = matrix_lakes_cor,  p.mat = matrix_lakes_p, method = "number",type = "upper",tl.cex = 1.5, tl.col = "black", number.cex = 1.5,hclust="ward",cl.pos = "n", sig.level = 0.99)
dev.off()


# lasso log RR -----

lasso_RRlog <- milasso(RRlogmice,M,"logRR",c("logRRn","logBdgT"))
write_xlsx(lasso_RRlog$lasso_coef,"8.version_control/lasso_coef_logRR.xlsx")

residuals <- lakes73log$logRR - lasso_RRlog$lasso_pred$pooled
x <- rnorm(73,0,1)

qqplot(x,residuals,xlab = "Normal distribution",ylab = "residuals",main = "Q-Q plot")
plot(lasso_RRlog$lasso_pred$pooled,residuals)+abline(h=0)
plot(ecdf(residuals))


plot(lakes73log$logRR,lasso_RRlog$lasso_pred$pooled)+abline(a=0,b=1)
rmse(lakes73log$logRR,lasso_RRlog$lasso_pred$pooled)

lm_RRlog <- with(RRlogmice,lm(logRR~logCN+Cells+c_CO2+SARuv+logFe+logDP+SUVA))
lm_RRlog_pooled <- pool(lm_RRlog) %>% summary()

lm_RRlog_predict <- with(RRlogmice,predict(lm(logRR~logCN+Cells+c_CO2+SARuv+logFe+logDP+SUVA)))
RRlog_predicted <- list.cbind(lm_RRlog_predict$analyses) %>% rowMeans()

rmse(lakes73log$logRR,RRlog_predicted)





# pearson log RR -----
predvar <- names(lakes73log)[-which(names(lakes73log) %in% c("logRR","logRRn"))]
pearson_pooled_RRlog <- predvar %>% as.data.frame() %>% setNames("param")
estimates <- c()
p_values <- c()
for(y in 1:length(predvar)){
  cor_formula <- as.formula(paste("~logRR+",predvar[y],sep=""))
  for(x in 1:M){
    df <- complete(RRlogmice,x) %>% dplyr::select(!c("logRRn"))
    dfcor <- cor.test(formula = cor_formula, data = df)
    estimates[x] <- dfcor$estimate[[1]]
    p_values[x] <- dfcor$p.value[[1]]
  }
  
  cor_pooled <- pool.scalar(estimates,p_values,n=73,k=1) 
  pearson_pooled_RRlog$r[y] <- cor_pooled$qbar
  pearson_pooled_RRlog$p_value[y] <- cor_pooled$ubar
  
}
pearson_pooled_RRlog <- pearson_pooled_RRlog[order(pearson_pooled_RRlog$r,decreasing = T),]
g <- ggplot(data = filter(pearson_pooled_RRlog,p_value <= 0.05),aes(x=reorder(param,order(r,decreasing=T)),y=r,label=param))+
  geom_col(aes(fill=p_value))+ 
  geom_text(angle=90,hjust=0,nudge_y = 0.01,size=7)+
  theme_bw(base_size=24)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.grid.minor = element_blank())+
  scale_fill_gradientn(colors = c("skyblue","dodgerblue4","firebrick"),
                       breaks = c(0,0.01,0.025, 0.05),labels = c("0","0.01","0.025","0.05"),limits = c(0,0.05), 
                       guide="legend")+
  labs(x = "", y= "r", title = "Pearson correlation coefficient with RR log") + ylim(-0.5,1)
print(g)
ggsave("8.version_control/pearson_cor_RRlog.png",width= 8,height = 5)


png("8.version_control/RRlog_cor.png",width=700,height = 700)
covar_RR <- filter(pearson_pooled_RRlog,p_value <= 0.05) %>% pull("param")
RR_cor <- micombine.cor(RRlogmice,variables = c("logRR",covar_RR), method="pearson") 
matrix_RR_cor <- attr(RR_cor,"r_matrix")
colnames(matrix_RR_cor) <- c("log(DN)","log(DP)", "log(SO4)","log(Na)","log(K)","log(Mg)","log(RR)","log(C:N)","log(C:P)")
rownames(matrix_RR_cor) <- c("log(DN)","log(DP)", "log(SO4)","log(Na)","log(K)","log(Mg)","log(RR)","log(C:N)","log(C:P)")
matrix_RR_p <-attr(RR_cor,"p_value")
corrplot::corrplot(corr = matrix_RR_cor,  p.mat = matrix_RR_p, method = "number",type = "upper",tl.cex = 1.5, tl.col = "black", number.cex = 1.5,hclust="ward",cl.pos = "n")
dev.off()

# Spearman log RR -----
spearman_pooled_logRR <- predvar %>% as.data.frame() %>% setNames("param")
estimates <- c()
p_values <- c()
for(y in 1:length(predvar)){
  cor_formula <- as.formula(paste("~logRR+",predvar[y],sep=""))
  for(x in 1:M){
    df <- complete(RRlogmice,x)
    dfcor <- cor.test(formula = cor_formula, data = df, method = "spearman",exact=F)
    estimates[x] <- dfcor$estimate[[1]]
    p_values[x] <- dfcor$p.value[[1]]
  }
  
  cor_pooled <- pool.scalar(estimates,p_values,n=73,k=1) 
  spearman_pooled_logRR$r[y] <- cor_pooled$qbar
  spearman_pooled_logRR$p_value[y] <- cor_pooled$ubar
  
}

spearman_pooled_logRR <- spearman_pooled_logRR[order(spearman_pooled_logRR$r,decreasing = T),]

g <- ggplot(data = filter(spearman_pooled_logRR,p_value <= 0.05),aes(x=reorder(param,order(r,decreasing=T)),y=r,label=reorder(param,order(r,decreasing=T))))+
  geom_col(aes(fill=p_value))+ 
  geom_text(angle=90,hjust=0,nudge_y = 0.01,size=7)+
  theme_bw(base_size=28)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.grid.minor = element_blank())+
  scale_fill_gradientn(colors = c("skyblue","dodgerblue4","firebrick"),
                       breaks = c(0,0.01,0.025, 0.05),labels = c("0","0.01","0.025","0.05"),limits = c(0,0.05), 
                       guide="legend")+
  xlab("")+ ylim(-1,1)
print(g)
ggsave("8.version_control/spearman_RRlog.png",width=8,height = 5)

png("8.version_control/logRR_scorplot.png",width = 600, height = 600)
covar_RR <- filter(spearman_pooled_logRR,p_value <= 0.05) %>% pull("param")
RR_scor <- micombine.cor(RRlogmice,variables = c("logRR",covar_RR), method="spearman") 
corrplot::corrplot(corr = attr(RR_scor,"r_matrix"), p.mat = attr(RR_scor,"p_matrix"),sig.level = 0.005, order="AOE",
                   method = "number",type = "upper",tl.cex = 1.5, tl.col = "black", number.cex = 1.5,cl.pos = "n")
dev.off()

# lasso log BdgT ----
lasso_BdgTlog <- milasso(RRlogmice,M,"logBdgT",c("logRR","logRRn"))
write_xlsx(lasso_BdgTlog$lasso_coef, "8.version_control/lasso_logBdgT.xlsx")

plot(lakes73log$BdgT,lasso_BdgTlog$lasso_pred$pooled)+abline(a=0,b=1)
rmse(lakes73log$logBdgT,lasso_BdgTlog$lasso_pred$pooled)

lm_BdgTlog <- with(RRlogmice,lm(logBdgT~logDOC+c_O2+pH+logAlkalinity+logDN+Cells+SARuv+logFe+logCO2))
lm_BdgTlog_pooled <- pool(lm_BdgTlog) %>% summary()

lm_BdgTlog_predict <- with(RRlogmice,predict(lm(logBdgT~logDOC+c_O2+pH+logAlkalinity+logDN+Cells+SARuv+logFe+logCO2)))
BdgTlog_predicted <- list.cbind(lm_BdgTlog_predict$analyses) %>% rowMeans()

rmse(lakes73log$logBdgT,BdgTlog_predicted)


# pearson log BdgT -----
predvar <- names(lakes73log)[-which(names(lakes73log) %in% c("logBdgT","logRRn"))]
pearson_pooled_BdgTlog <- predvar %>% as.data.frame() %>% setNames("param")
estimates <- c()
p_values <- c()
for(y in 1:length(predvar)){
  cor_formula <- as.formula(paste("~logBdgT+",predvar[y],sep=""))
  for(x in 1:M){
    df <- complete(RRlogmice,x) 
    dfcor <- cor.test(formula = cor_formula, data = df)
    estimates[x] <- dfcor$estimate[[1]]
    p_values[x] <- dfcor$p.value[[1]]
  }
  
  cor_pooled <- pool.scalar(estimates,p_values,n=73,k=1) 
  pearson_pooled_BdgTlog$r[y] <- cor_pooled$qbar
  pearson_pooled_BdgTlog$p_value[y] <- cor_pooled$ubar
  
}
pearson_pooled_BdgTlog <- pearson_pooled_BdgTlog[order(pearson_pooled_BdgTlog$r,decreasing = T),]
g <- ggplot(data = filter(pearson_pooled_BdgTlog,p_value <= 0.05),aes(x=reorder(param,order(r,decreasing=T)),y=r,label=param))+
  geom_col(aes(fill=p_value))+ 
  geom_text(angle=90,hjust=0,nudge_y = 0.01,size=7)+
  theme_bw(base_size=24)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.grid.minor = element_blank())+
  scale_fill_gradientn(colors = c("skyblue","dodgerblue4","firebrick"),
                       breaks = c(0,0.01,0.025, 0.05),labels = c("0","0.01","0.025","0.05"),limits = c(0,0.05), 
                       guide="legend")+
  labs(x = "", y= "r", title = "Pearson correlation coefficient with BdgT log") + ylim(-0.5,1)
print(g)
ggsave("8.version_control/pearson_cor_BdgTlog.png",width= 8,height = 5)


png("8.version_control/BdgTlog_cor.png",width=700,height = 700)
covar_BdgT <- filter(pearson_pooled_BdgTlog,p_value <= 0.05) %>% pull("param")
BdgT_cor <- micombine.cor(RRlogmice,variables = c("logBdgT",covar_BdgT), method="pearson")
matrix_BdgT_cor <- attr(BdgT_cor,"r_matrix")
colnames(matrix_BdgT_cor) <- c("log(DN)","log(Alkalinity)","log(SO4)","log(K)","log(Ca)","log(Mg)","log(O2)","log(BdgT)","pH","c_O2")
rownames(matrix_BdgT_cor) <- c("log(DN)","log(Alkalinity)","log(SO4)","log(K)","log(Ca)","log(Mg)","log(O2)","log(BdgT)","pH","c_O2")
matrix_BdgT_p <-attr(BdgT_cor,"p_value")
corrplot::corrplot(corr = matrix_BdgT_cor,  p.mat = matrix_BdgT_p, method = "number",type = "upper",tl.cex = 1.5, tl.col = "black", number.cex = 1.5,hclust="ward",cl.pos = "n")
dev.off()

# Spearman log BdgT -----
spearman_pooled_logBdgT <- predvar %>% as.data.frame() %>% setNames("param")
estimates <- c()
p_values <- c()
for(y in 1:length(predvar)){
  cor_formula <- as.formula(paste("~logBdgT+",predvar[y],sep=""))
  for(x in 1:M){
    df <- complete(RRlogmice,x)
    dfcor <- cor.test(formula = cor_formula, data = df, method = "spearman",exact=F)
    estimates[x] <- dfcor$estimate[[1]]
    p_values[x] <- dfcor$p.value[[1]]
  }
  
  cor_pooled <- pool.scalar(estimates,p_values,n=73,k=1) 
  spearman_pooled_logBdgT$r[y] <- cor_pooled$qbar
  spearman_pooled_logBdgT$p_value[y] <- cor_pooled$ubar
  
}

spearman_pooled_logBdgT <- spearman_pooled_logBdgT[order(spearman_pooled_logBdgT$r,decreasing = T),]

g <- ggplot(data = filter(spearman_pooled_logBdgT,p_value <= 0.05),aes(x=reorder(param,order(r,decreasing=T)),y=r,label=reorder(param,order(r,decreasing=T))))+
  geom_col(aes(fill=p_value))+ 
  geom_text(angle=90,hjust=0,nudge_y = 0.01,size=7)+
  theme_bw(base_size=28)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.grid.minor = element_blank())+
  scale_fill_gradientn(colors = c("skyblue","dodgerblue4","firebrick"),
                       breaks = c(0,0.01,0.025, 0.05),labels = c("0","0.01","0.025","0.05"),limits = c(0,0.05), 
                       guide="legend")+
  xlab("")+ ylim(-1,1)
print(g)
ggsave("8.version_control/spearman_BdgTlog.png",width=8,height = 5)

png("8.version_control/logBdgT_scorplot.png",width = 600, height = 600)
covar_BdgT <- filter(spearman_pooled_logBdgT,p_value <= 0.05)%>% filter(abs(r) > 0.3) %>% pull("param")
BdgT_scor <- micombine.cor(RRlogmice,variables = c("logBdgT",covar_BdgT), method="spearman") 
corrplot::corrplot(corr = attr(BdgT_scor,"r_matrix"), p.mat = attr(BdgT_scor,"p_matrix"),sig.level = 0.005, order="AOE",
                   method = "number",type = "upper",tl.cex = 2, tl.col = "black", number.cex = 2,cl.pos = "n")
dev.off()

# lasso log RRn -----
lasso_RRnlog <- milasso(RRlogmice,M,"logRRn","logRR")
write_xlsx(lasso_RRnlog$lasso_coef,"8.version_control/lasso_coef_logRRn.xlsx")

# linear regression log RRn -----
plot(lakes73log$logRRn,lasso_RRnlog$lasso_pred$pooled)+abline(a=0,b=1)
rmse(lakes73log$logBdgT,lasso_RRnlog$lasso_pred$pooled)

lm_RRnlog <- with(RRlogmice,lm(logRRn~SR+c_O2+pH+c_CO2+SARuv+logMg+logFe+logNa+logCH4+logDP+logDN+SUVA))
lm_RRnlog_pooled <- pool(lm_RRnlog) %>% summary()

lm_RRnlog_predict <- with(RRlogmice,predict(lm(logRRn~SR+c_O2+pH+c_CO2+SARuv+logMg+logFe+logNa+logCH4+logDP+logDN+SUVA)))
RRnlog_predicted <- list.cbind(lm_RRnlog_predict$analyses) %>% rowMeans()

plot(lakes73log$logRRn,RRnlog_predicted)+abline(a=0,b=1)
rmse(lakes73log$logRRn,RRnlog_predicted)

# Pearson log RRn -----
predvar <- names(lakes73log)[-which(names(lakes73log) == "logRRn")]
pearson_pooled_RRnlog <- predvar %>% as.data.frame() %>% setNames("param")
estimates <- c()
p_values <- c()
for(y in 1:length(predvar)){
  cor_formula <- as.formula(paste("~logRRn+",predvar[y],sep=""))
  for(x in 1:M){
    df <- complete(RRlogmice,x)
    dfcor <- cor.test(formula = cor_formula, data = df)
    estimates[x] <- dfcor$estimate[[1]]
    p_values[x] <- dfcor$p.value[[1]]
  }
  
  cor_pooled <- pool.scalar(estimates,p_values,n=73,k=1) 
  pearson_pooled_RRnlog$r[y] <- cor_pooled$qbar
  pearson_pooled_RRnlog$p_value[y] <- cor_pooled$ubar
  
}
pearson_pooled_RRnlog <- pearson_pooled_RRnlog[order(pearson_pooled_RRnlog$r,decreasing = T),]
g <- ggplot(data = filter(pearson_pooled_RRnlog,p_value <= 0.05),aes(x=reorder(param,order(r,decreasing=T)),y=r,label=param))+
  geom_col(aes(fill=p_value))+ 
  geom_text(angle=90,hjust=0,nudge_y = 0.01,size=7)+
  theme_bw(base_size=24)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.grid.minor = element_blank())+
  scale_fill_gradientn(colors = c("skyblue","dodgerblue4","firebrick"),
                       breaks = c(0,0.01,0.025, 0.05),labels = c("0","0.01","0.025","0.05"),limits = c(0,0.05), 
                       guide="legend")+
  labs(x = "", y= "r", title = "Pearson correlation coefficient with RR log") + ylim(-1,1)
print(g)
ggsave("8.version_control/pearson_RRnlog.png",width= 8,height = 5)


png("8.version_control/RRnlog_cor.png",width=700,height = 700)
covar_RRn <- filter(pearson_pooled_RRnlog,p_value <= 0.05) %>% filter(abs(r) > 0.4) %>% pull("param")
RRn_cor <- micombine.cor(RRlogmice,variables = c("logRRn",covar_RRn), method="pearson") 
matrix_RRn_cor <- attr(RRn_cor,"r_matrix")
colnames(matrix_RRn_cor) <- c("log(DOC)","log(DN)","log(EC)","log(Cl)","log(B)","log(Na)","log(K)","log(Mg)","log(CH4)","log(RR)", "log(RRn)","SUVA","SR")
rownames(matrix_RRn_cor) <- c("log(DOC)","log(DN)","log(EC)","log(Cl)","log(B)","log(Na)","log(K)","log(Mg)","log(CH4)","log(RR)", "log(RRn)","SUVA","SR")
matrix_RRn_p <-attr(RRn_cor,"p_value")
corrplot::corrplot(corr = matrix_RRn_cor,  p.mat = matrix_RRn_p, method = "number",type = "upper",tl.cex = 1.3, tl.col = "black", number.cex = 1.3,hclust="ward",cl.pos = "n")
dev.off()
