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
library(gridExtra)
library(grid)


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
  
  lasso_raw <- mice.object$data[,resp.var]
  
  n <- dim(mice.object$data)[1] 
  lasso_coef <- c("Intercept",names(dplyr::select(mice.object$data,!c(all_of(resp.var),all_of(excluded.var))))) %>% as.data.frame() %>% setNames("param")
  
  lasso_pred <- c(rep(NA,n)) %>% as.data.frame()
  lasso_res <- c(rep(NA,n)) %>% as.data.frame()
  
  lasso_lambda <- matrix(data=NA,nrow = 100, ncol=M) %>% as.data.frame()
  lasso_cvm <- matrix(data=NA,nrow = 100, ncol=M) %>% as.data.frame()
  lasso_cvlo <- matrix(data=NA,nrow = 100, ncol=M) %>% as.data.frame()
  lasso_cvup <- matrix(data=NA,nrow = 100, ncol=M) %>% as.data.frame()
  lambda_min <- c()
  lambda_1se <- c()
  
  
  for (z in 1:M){
    lasso_coef %>% tibble::add_column(z=NA)
    data <- complete(mice.object,z) %>% dplyr::select(!excluded.var)
    set.seed(5)
    lasso.formula <- as.formula(paste(resp.var,"~.",sep=""))
    fit <- cv.glmnet(formula = lasso.formula,data = data)
    coef_fit <- fit %>% coef(s="lambda.min") %>% summary()
    l <- length(fit$lambda)
    
    for (w in 1:dim(coef_fit)[1]){
      lasso_coef[coef_fit[w,1],z+1] <- coef_fit[w,3]
    }
    
    newmat <- complete(mice.object,z) %>% dplyr::select(!c(resp.var,excluded.var)) 
    lasso_pred[,z] <- predict(fit,newmat[1:n,],s="lambda.min")
    lasso_res[,z] <- data[,resp.var] - lasso_pred[,z]
    lasso_lambda[1:l,z] <- fit$lambda
    lasso_cvm[1:l,z] <- fit$cvm
    lasso_cvlo[1:l,z] <- fit$cvlo
    lasso_cvup[1:l,z] <- fit$cvup
    lambda_min[z] <- fit$lambda.min
    lambda_1se[z] <- fit$lambda.1se
  }
  
  # pools lasso_pred and lasso_res
  lasso_pred$pooled <- rowMeans(lasso_pred)
  lasso_res$pooled <- rowMeans(lasso_res)

  # pools lambdas
  
  lasso_lambda <- lasso_lambda[rowSums(is.na(lasso_lambda)) != ncol(lasso_lambda), ] # removes rows with only NA in lasso_lambda
  k <- dim(lasso_lambda)[1]
  
  lasso_cv <- c(rep(NA,k)) %>% as.data.frame()

  lasso_cv$lambda <- rowMeans(lasso_lambda,na.rm = T)
  
  lasso_cvm <- lasso_cvm[rowSums(is.na(lasso_cvm)) != ncol(lasso_cvm), ] # removes rows with only NA in lasso_cvm
  lasso_cv$cvm <- rowMeans(lasso_cvm,na.rm = T)
  
  lasso_cvlo <- lasso_cvlo[rowSums(is.na(lasso_cvlo)) != ncol(lasso_cvlo), ] # removes rows with only NA in lasso_cvlo
  lasso_cv$cvlo <- rowMeans(lasso_cvlo,na.rm = T)
  
  lasso_cvup <- lasso_cvup[rowSums(is.na(lasso_cvup)) != ncol(lasso_cvup), ] # removes rows with only NA in lasso_cvup
  lasso_cv$cvup <- rowMeans(lasso_cvup,na.rm = T)
  
  # prepares lasso_coef for final file
  lasso_coef$pooled <- NA
  lasso_coef$number <- NA
  
    for(x in 1:length(lasso_coef$param)){
    estimates <- lasso_coef[x,] %>% dplyr::select(c(1:M+1)) %>% as.numeric()
    lasso_coef$pooled[x] <- mean(estimates,na.rm = T)
    lasso_coef$number[x] <- M - sum(is.na(estimates))
  }

  final <- list(lasso_raw, lasso_coef,lasso_pred,lasso_res,lasso_lambda,lasso_cv,lambda_min,lambda_1se)
  names(final) <- c("lasso_raw","lasso_coef","lasso_pred","lasso_res","lasso_lambda","lasso_cv","lambda_min","lambda_1se")
  return(final)

  
}

lasso_res_plot <- function(lasso_list, title = "residual plot"){
  
  w <- lasso_list$lasso_raw
  x <- lasso_list$lasso_pred$pooled
  j <- order(x)
  y <- lasso_list$lasso_res$pooled
  
  # plots cross-validation
  g1 <- ggplot(lasso_list$lasso_cv,aes(x=log(lambda)))+geom_point(aes(y=cvm),col="red")+
    geom_errorbar(aes(ymin=cvlo,ymax=cvup),col="gray")+
    geom_vline(xintercept=log(mean(lasso_list$lambda_min)),linetype = "dashed")+
    geom_label(aes(x=log(mean(lasso_list$lambda_min)),y=max(cvup)/2,label="lambda min",angle=90))+
    geom_vline(xintercept=log(mean(lasso_list$lambda_1se)), linetype = "dashed")+
    geom_label(aes(x=log(mean(lasso_list$lambda_1se)),y=max(cvup)/3,label="lambda 1se",angle=90))+
    labs(x="log(lambda)",y="Mean squared error",title = "Cross-validated lambda")+
    theme_bw(base_size=15)+ theme(panel.grid = element_line(color="gray95"))
  
  # Fitted vs observed
  g2 <- qplot(w,x)+
    labs(x = "Observed values", y = "Fitted values" , title = "Fitted values vs. observations")+
    annotate("label",x = max(w) - (max(w)-min(w))/6, y = min(x) + (max (x)-min(x))/6, label = paste("MAE",round(mae(w,x),2),sep = " = "), size = 5)+
    theme_bw(base_size=15)+ theme(panel.grid = element_line(color="gray95"))
  
  # QQ plot
  g3 <- ggplot(lasso_list$lasso_res,aes(sample=pooled)) + stat_qq() + stat_qq_line(col = "gray") +
    labs (x = "Theoretical quantiles", y = "Sample quantiles" , title = "Normal Q-Q plot")+
    theme_bw(base_size = 15) +  theme(panel.grid = element_line(color="gray95"))   
  
  # residual vs predicted
  
  z <- loess(y~x)
  
  g4 <- qplot(x,y,xlab = "Predicted", ylab = "Residuals") + geom_hline (yintercept = 0,col = "gray") + 
    geom_smooth(method = "loess", se = F, col = "red")+
#    geom_line(aes(x=x[j],y = z$fitted[j]),col="red",lwd=1)+
    labs (x = "Fitted values", title = "Residuals vs fitted")+
    theme_bw(base_size = 15) + theme(panel.grid = element_line(color="gray95"))
  
  grid.arrange(g2,g1,g3,g4,ncol = 2, top = textGrob(title,gp = gpar(fontsize = 25, font = 3),vjust = 0.5))
}

lasso_simple_plot <- function(lasso_list, title = ""){
  
  w <- lasso_list$lasso_raw
  x <- lasso_list$lasso_pred$pooled
  j <- order(x)
  y <- lasso_list$lasso_res$pooled

  # Fitted vs observed
  g1 <- qplot(w,x,col=y)+
    labs(x = "Observed values", y = "Fitted values" , title = paste("Fitted values vs. observations",title, sep = " "),col="Residuals")+
    geom_abline(intercept = 0, slope = 1, col = "gray")+
    scale_color_viridis_c(direction = -1, end = 1 ,begin = 0)+
    annotate("label",x = max(w) - (max(w)-min(w))/6, y = min(x) + (max (x)-min(x))/6, label = paste("MAE",round(mae(w,x),2),sep = " = "), size = 5)+
    theme_bw(base_size=15)+ theme(panel.grid = element_line(color="gray95"))
  plot(g1)

}



gauss_lasso <- function (mice.object,resp.var,lasso_list){
  
  title <- paste("Residual plots linear model", resp.var, sep = " ")
  
  lm_obs <- mice.object$data[,resp.var]
  n <- length(lm_obs)
  
  lasso_param <- lasso_list$lasso_coef$param %>% as.data.frame %>% setNames("lasso_param")
  covar_index <- which(lasso_list$lasso_coef$number > mice.object$m/2)
  covar_lm <- lasso_list$lasso_coef$param[covar_index]
  p <- length(covar_lm)
  
  lm_fm <- paste(resp.var,paste(covar_lm[-1],collapse = "+"), sep = "~")
  
  lm_df <- with(mice.object,lm(as.formula(lm_fm))) %>% pool() %>% summary() %>% as.data.frame() %>% select(c("term","estimate","p.value"))
  levels(lm_df$term)[which(levels(lm_df$term) == "(Intercept)")] <- "Intercept"
  fit_print <- merge(lasso_param,lm_df, by.x= "lasso_param",by.y = "term", all.x = T) %>% mutate(Response = resp.var)
  assign(paste("lm_fit_",resp.var,sep=""), fit_print,envir = .GlobalEnv)
  write_xlsx(fit_print,paste("8.version_control/gauss_lasso_coef_",resp.var,".xlsx",sep=""))
  
  lm_pred_with <- with(mice.object,predict(lm(as.formula(lm_fm))))
  lm_pred <- list.cbind(lm_pred_with$analyses) %>% rowMeans()
  
  lm_res_with <- with(mice.object,residuals(lm(as.formula(lm_fm))))
  lm_res <- list.cbind(lm_res_with$analyses) %>% rowMeans()
  
  lm_rstandard_with <- with(mice.object,rstandard(lm(as.formula(lm_fm))))
  lm_rstandard <- list.cbind(lm_rstandard_with$analyses) %>% rowMeans()
  
  lm_lev_with <- with(mice.object,hatvalues(lm(as.formula(lm_fm))))
  lm_lev <- list.cbind(lm_lev_with$analyses) %>% rowMeans()
  
  lm_cook_with <- with(mice.object,cooks.distance(lm(as.formula(lm_fm))))
  lm_cook <- list.cbind(lm_cook_with$analyses) %>% rowMeans()
  D <- 4/(n-p-1)
  outliers <- which(lm_cook > D)
  names(outliers) <- lakes73$Lake_name[outliers]
  
  
  # fitted vs observed
  g1 <- qplot(lm_obs,lm_pred)+
    labs(x = "Observed values", y = "Fitted values" , title = "Fitted values vs. observations")+
    annotate("label",x = max(lm_obs) - (max(lm_obs)-min(lm_obs))/6, y = min(lm_pred) + (max (lm_pred)-min(lm_obs))/6, label = paste("MAE",round(mae(lm_obs,lm_pred),2),sep = " = "), size = 5)+
    theme_bw(base_size=15)+ theme(panel.grid = element_line(color="gray95"))
  
  #Residuals vs fitted
  g2 <- qplot(lm_pred,lm_res)+geom_smooth(method = "loess", size = 1, se = F, col = "red")+
    labs(x="Fitted", y = "Residuals", title = "Residuals vs fitted")+
    theme_bw(base_size = 15)+theme(panel.grid = element_line(color = "gray95"))
  
  # standardised residuals 
  g3 <- qplot(lm_pred,sqrt(abs(lm_rstandard)))+geom_smooth(method = "loess", size = 1, se = F, col = "red")+
    labs(x="Fitted", y = "sqrt Standardised residuals", title = "Scale-location")+
    theme_bw(base_size = 15)+theme(panel.grid = element_line(color = "gray95"))
  
  #QQplot
  g4 <- qplot(sample=lm_res)+stat_qq_line(col = "gray")+
    labs (x = "Theoretical quantiles", y = "Sample quantiles" , title = "Normal Q-Q plot")+
    theme_bw(base_size = 15) +  theme(panel.grid = element_line(color="gray95"))  
  
  # Residuals vs leverage
  g5 <- qplot(lm_lev,lm_rstandard)+geom_smooth(method = "loess", size = 1, se = F, col = "red")+
    labs(x="Leverage", y = "Standardised residuals", title = "Residuals vs leverage")+
    theme_bw(base_size = 15)+theme(panel.grid = element_line(color = "gray95"))+ 
    geom_line(aes(x= lm_lev,y = sqrt(D*p*(1-lm_lev)/lm_lev)),col="darkblue",lty = 2)+
    geom_line(aes(x = lm_lev, y = -sqrt(D*p*(1-lm_lev)/lm_lev)),col="darkblue",lty=2)+
    geom_label_repel(aes(x=lm_lev[outliers],y=lm_rstandard[outliers]),label=names(outliers),nudge_y = 0.5,nudge_x = 0.1, size = 3.5)+
    ylim(min(lm_rstandard)-sd(lm_rstandard),max(lm_rstandard)+sd(lm_rstandard))+
    annotate("label",x = max(lm_lev)*0.8 ,y= min(lm_rstandard),label="Cook's distance (n-p-1)",col="darkblue",size = 3.5)
  
  g6 <- qplot(lm_pred,lm_rstandard)+geom_smooth(method = "loess", size = 1, se = F, col = "red")+
    labs(x="Fitted", y = "Residuals", title = "Residuals vs fitted")+
    theme_bw(base_size = 15)+theme(panel.grid = element_line(color = "gray95"))
  
  grid.arrange(g1,g2,g4,g5,ncol = 2, top = textGrob(title,gp = gpar(fontsize = 25, font = 3)))
  
  g7 <- qplot(lm_obs,lm_pred,col=lm_rstandard)+
    labs(x = "Observed values", y = "Fitted values" , title = paste("Fitted values vs. observations",resp.var,sep = " "),col="Standardised Residuals")+
    geom_abline(intercept = 0, slope = 1, col = "gray")+
    scale_color_viridis_c(direction = -1, end = 1 ,begin = 0)+
    annotate("label",x = max(w) - (max(w)-min(w))/6, y = min(x) + (max (x)-min(x))/6, label = paste("MAE",round(mae(w,x),2),sep = " = "), size = 5)+
    theme_bw(base_size=20)+ theme(panel.grid = element_line(color="gray95"))
  plot(g7)
  
}  



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
logcovariates <- c("logRR","logRRn","logBdgT","logDOC","logCN","logDP","logEC","logFe","pH","Cells","SUVA","SARuv","logO2","logN2O","p_O2","p_CO2","c_CO2","p_CH4")
RRlogmice <- select(lakes73log,logcovariates) %>% mice(method = "cart", m = M)

lasso_var <- c("Intercept","logDOC","logDP","logEC","logFe","logO2","logN2O","logCN","pH","Cells","SUVA","SARuv","p_O2","p_CO2","c_CO2","p_CH4")

# correlogram for RRlogmice -----
png("8.version_control/correlogram.png", width= 1000, height = 1000)
lakes_cor <- micombine.cor(RRlogmice, method="pearson") 
matrix_lakes_cor <- attr(lakes_cor,"r_matrix")
colnames(matrix_lakes_cor) <- c("log(RR)","log(RRn)", "log(BdgT)","log(DOC)","log(C:N)","log(DP)","log(EC)","log(Fe)","pH","Cells","SUVA","SARuv","log(O2)","log(N2O)","p_O2","p_CO2","c_CO2","p_CH4")
rownames(matrix_lakes_cor) <- c("log(RR)","log(RRn)", "log(BdgT)","log(DOC)","log(C:N)","log(DP)","log(EC)","log(Fe)","pH","Cells","SUVA","SARuv","log(O2)","log(N2O)","p_O2","p_CO2","c_CO2","p_CH4")
matrix_lakes_p <-attr(lakes_cor,"p_value")
corrplot::corrplot(corr = matrix_lakes_cor,  p.mat = matrix_lakes_p, method = "number",type = "upper",tl.cex = 1.5, tl.col = "black", number.cex = 1.5,hclust="single",cl.pos = "n", sig.level = 0.99)
dev.off()


# lasso  -----

lasso_RRlog <- milasso(RRlogmice,M,"logRR",c("logRRn","logBdgT"))
write_xlsx(lasso_RRlog$lasso_coef,"8.version_control/lasso_coef_logRR.xlsx")

lasso_BdgTlog <- milasso(RRlogmice,M,"logBdgT",c("logRR","logRRn"))
write_xlsx(lasso_BdgTlog$lasso_coef, "8.version_control/lasso_logBdgT.xlsx")

lasso_RRnlog <- milasso(RRlogmice,M,"logRRn",c("logRR","logBdgT"))
write_xlsx(lasso_RRnlog$lasso_coef,"8.version_control/lasso_coef_logRRn.xlsx")

df1 <- select(lasso_RRlog$lasso_coef,c("param","pooled","number")) %>% mutate(Response = "logRR")
df2 <- select(lasso_BdgTlog$lasso_coef,c("param","pooled","number")) %>% mutate(Response = "logBdgT")
df3 <- select(lasso_RRnlog$lasso_coef,c("param","pooled","number")) %>% mutate(Response = "logRRn")
long_lasso <- rbind(df1,df2,df3) %>% filter(number > M/2) %>% filter(param != "Intercept")


# plots lasso on pdf -----
pdf("8.version_control/lasso_res_plot.pdf",height = 10, width = 15)
lasso_res_plot(lasso_RRlog, title = "Residual plots lasso regression RR")
lasso_res_plot(lasso_BdgTlog, title = "Residual plots lasso regression BdgT")
lasso_res_plot(lasso_RRnlog, title = "Residual plots lasso regression RRn")

dev.off()

pdf("8.version_control/lasso_model_plot.pdf",height = 8,width = 10)

lasso_simple_plot(lasso_RRlog,title = "log(RR)")
lasso_simple_plot(lasso_RRnlog,title = "log(RRn)")
lasso_simple_plot(lasso_BdgTlog, title = "log(BdgT)")

ggplot(long_lasso)+geom_col(aes(x=param,y=pooled,fill=Response))+facet_grid(rows = vars(Response),scales = "free_y")+
  labs(x="",y="Pooled coefficient")+
  theme_light(base_size = 15)+
  theme(legend.position = "none",axis.text.x = element_text(angle = 45,hjust = 1),panel.grid.minor = element_line(colour = NA))+
  scale_fill_viridis_d(end = 0.8)

ggplot(long_lasso)+geom_col(aes(x=param,y=pooled,fill=Response))+facet_grid(rows = vars(Response))+
  labs(x="",y="Pooled coefficient")+
  theme_light(base_size = 15)+
  theme(legend.position = "none",axis.text.x = element_text(angle = 45,hjust = 1),panel.grid.minor = element_line(colour = NA))+
  scale_fill_viridis_d(end = 0.8)

dev.off()

  

# Linear model  ----- 

pdf("8.version_control/gauss_lasso.pdf",height = 10,width = 15)
gauss_lasso(RRlogmice,"logRR",lasso_RRlog)
gauss_lasso(RRlogmice,"logBdgT",lasso_BdgTlog)
gauss_lasso(RRlogmice,"logRRn",lasso_RRnlog)

dev.off()

pdf("8.version_control/gauss_lasso_model.pdf",height = 8, width = 10)
long_lm <- rbind(lm_fit_logRR,lm_fit_logRRn,lm_fit_logBdgT) %>% filter(lasso_param != "Intercept") %>% filter(p.value < 0.05)

ggplot(long_lm)+geom_col(aes(x=lasso_param,y=estimate,fill=Response))+facet_grid(rows = vars(Response),scales = "free_y")+
  labs(x="",y="Pooled coefficient")+
  theme_light(base_size = 15)+
  theme(legend.position = "none",axis.text.x = element_text(angle = 45,hjust = 1),panel.grid.minor = element_line(colour = NA))+
  scale_fill_viridis_d(end = 0.8)

ggplot(long_lm)+geom_col(aes(x=lasso_param,y=estimate,fill=Response))+facet_grid(rows = vars(Response))+
  labs(x="",y="Pooled coefficient")+
  theme_light(base_size = 15)+
  theme(legend.position = "none",axis.text.x = element_text(angle = 45,hjust = 1),panel.grid.minor = element_line(colour = NA))+
  scale_fill_viridis_d(end = 0.8)

dev.off()

# removing one outlier -----
set.seed(5)
RRlogmice2 <- select(lakes73log[-outliers,],logcovariates) %>%  mice(method = "cart", m = M)

lasso_RRlog2 <- milasso(RRlogmice2,M,"logRR",c("logRRn","logBdgT"))
write_xlsx(lasso_RRlog2$lasso_coef,"8.version_control/lasso_coef_logRR2.xlsx")

lasso_BdgTlog2 <- milasso(RRlogmice2,M,"logBdgT",c("logRR","logRRn"))
write_xlsx(lasso_BdgTlog2$lasso_coef, "8.version_control/lasso_logBdgT2.xlsx")

lasso_RRnlog2 <- milasso(RRlogmice2,M,"logRRn","logRR")
write_xlsx(lasso_RRnlog2$lasso_coef,"8.version_control/lasso_coef_logRRn2.xlsx")

pdf("8.version_control/lasso_res_plot2.pdf",height = 10, width = 15)
lasso_res_plot(lasso_RRlog2, title = "Residual plots lasso regression RR")
lasso_res_plot(lasso_BdgTlog2, title = "Residual plots lasso regression BdgT")
lasso_res_plot(lasso_RRnlog2, title = "Residual plots lasso regression RRn")

dev.off()


pdf("8.version_control/gauss_lasso2.pdf",height = 10,width = 15)
gauss_lasso(RRlogmice2,"logRR",lasso_RRlog2)
gauss_lasso(RRlogmice2,"logBdgT",lasso_BdgTlog2)
gauss_lasso(RRlogmice2,"logRRn",lasso_RRnlog2)
dev.off()
