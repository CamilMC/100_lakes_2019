# spectral slopes -----

getslope <- function(df,lambda1,lambda2){
  a1 <- df %>% filter(WL.nm == lambda1)
  a2 <- df %>% filter(WL.nm == lambda2)
  a <- rbind(a1,a2) %>% t() %>% as.data.frame() %>% tibble::rownames_to_column() %>% filter(rowname != "WL.nm") 
  a <- apply(a,2,as.numeric) %>% as.data.frame()
  a$s <- (-1/(lambda2-lambda1))*log(a[,3]/a[,2])
  a <- setNames(a,c("Lake_ID",paste("a",lambda1,sep=""),paste("a",lambda2,sep=""),paste("s",lambda1,lambda2,sep = "_")))
  return(a)
}

# histograms -----
# returns histograms of each of the columns of a dataframe

hist.df <- function(df){
  for (i in names(df)){
    print(ggplot(df)+ geom_histogram(aes(df[,i]),col="black", bins = 100,fill="white")+
            xlab(i)+theme_light(base_size=20))
  }
}

# boxplots -----
# returns boxplot of each of the columns of a dataframe + a list of the outliers in a list object called "out"

boxplots <- function(df,coef_out){
  library("ggrepel")
  
  out <- list()
  
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
  
  out.df <- data.frame(data = df$Sample_ID) %>% setNames("Sample_ID")
  out.df$outlier <- 0
  
  
  for(j in 1:length(out.df$Sample_ID)){
    param.out <- grep(out.df$Sample_ID[j],out)
    out.df$outlier[j] <- length(param.out) %>% as.numeric()
    out.df$param[j] <- names(out)[param.out] %>% toString()
  }
  
  
  out.df.param <- data.frame(data = names(out)) %>% setNames("param")
  out.df.param$outlier <- NA
  
  for(k in 1:length(out.df.param$param)){
    out.df.param$outlier[k] <- out[[k]]$Sample_ID %>% length()
    out.df.param$Sample_ID[k] <- out[[k]]$Sample_ID %>% toString()
  }
  
  out <<- out
  out.df <<- out.df
  out.df.param <<- out.df.param
  
}

# plot.cor -----
# plots correlations for a specific parameter in a dataframe
# returns a plot

plot.cor <- function(df,param,lim=c(-1,1)){
  df2 <- df
  name_df <- deparse(substitute(df))
  cor.matrix <- cor(df2,use="pairwise")
  cor.param <- cor.matrix[,which(colnames(cor.matrix)==param)] %>% sort(decreasing = T) %>% 
    as.data.frame() %>% tibble::rownames_to_column() %>% setNames(c("predictor","corr.param"))
  cor.param.p <- cor.param[-1,]
  g <- ggplot(cor.param.p,aes(x=reorder(predictor,order(cor.param.p$corr.param,decreasing=T)),y=corr.param,label=predictor))+
    geom_col(aes(fill=corr.param))+ 
    geom_text(angle=90,hjust=0,nudge_y = 0.01,size=6)+
    theme_bw(base_size=20)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    scale_fill_gradient(high="red",low="skyblue")+ylim(lim)+labs(x="",y=paste(name_df,param,sep="."),fill="correlation")
  #scale_fill_gradientn(colors = c("skyblue","dodgerblue4","firebrick"),values = c(0,0.01,0.05,1),
   #                    breaks = c(0,0.01,0.05,1),labels = c("0","0.01","0.05","1"),limits = c(0,1), guide="legend")+
    print(g)
  ggsave(plot = g, filename = paste("cor",name_df,param,"png",sep="."),device ="png",width=15,height=10)
}

# table.cor -----
# writes a table with the correlations between one specific parameter in one dataframe
# returns a table called "cor.param" 
# returns a csv file in the "5. 100 lakes" folder with the name of the df and of the parameter

table.cor <- function (df,param){
  #save name of the dataframe as a character
  df2 <- df
  name_df <- deparse(substitute(df))
  
  # correlation table with p values
  cor.param <-data.frame(predictor = rep(NA,length(names(df2)[-1])))
  cor.param$predictor <- names(df2)[-1] %>% as.vector()
  cor.param$corr <- sapply(cor.param[,1], function(x) cor.test(waterchem[[param]],waterchem[[x]])$estimate)
  cor.param$pvalue <- sapply(cor.param[,1], function(x) cor.test(waterchem[[param]],waterchem[[x]])$p.value)
  cor.param <- cor.param[order(cor.param$pvalue),]
  write.csv(cor.param,paste("5. 100 lakes/",name_df,".",param,".csv",sep = ""))
}


# print.cor -----
# prints both the correlations between one parameter and the rest of the df
# returns plots
# adds an element to the list "all.correlations"

print.cor <- function(df,param,title = T,cor_list = T){
  
  
  #save name of the dataframe as a character
  df2 <- df
  name_df <- paste(deparse(substitute(df)),param,sep=".")
  
  # correlation table with p values
  cor.param <-data.frame(predictor = rep(NA,length(names(df2)[-1])))
  cor.param$predictor <- names(df2)[-1] %>% as.vector()
  cor.param$corr <- sapply(cor.param[,1], function(x) cor.test(df2[[param]],df2[[x]])$estimate)
  cor.param$pvalue <- sapply(cor.param[,1], function(x) cor.test(df2[[param]],df2[[x]])$p.value)
  cor.param <- cor.param[order(cor.param$corr, decreasing = T),]
  
  #plot cor
  g <- ggplot(data = cor.param[-1,],aes(x=reorder(predictor,order(cor.param[-1,]$corr,decreasing=T)),
                                        y=corr,label=predictor))+
    geom_col(aes(fill=pvalue))+ 
    geom_text(angle=90,hjust=0,nudge_y = 0.01,size=6)+
    theme_bw(base_size=20)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    scale_fill_gradientn(colors = c("skyblue","dodgerblue4","firebrick"),values = c(0,0.01,0.05,1),
                         breaks = c(0,0.01,0.05,1),labels = c("0","0.01","0.05","1"),limits = c(0,1), guide="legend")+
    ylim(-1.1,1.1)
  if(title == T){
    g <- g + labs(x="",y=name_df,fill="p-value")
    print(g)
  } else if(title == F) {
    g <- g + labs(x="",y=paste("correlation",param,sep="_"),fill="p-value")
    print(g)
  }
  
  #print table
  if(cor_list == T){
    all.correlations[[param]] <<- cor.param
  } 
}

# print.cor.signif -----
# only prints the significant correlations on the graph
# returns graphs but does not change the list "all.correlations"

print.cor.signif <- function(df,param,title = T){
  
  #save name of the dataframe as a character
  df2 <- df
  name_df <- paste(deparse(substitute(df)),param,sep=".")
  
  # correlation table with p values
  cor.param <-data.frame(predictor = rep(NA,length(names(df2))))
  cor.param$predictor <- names(df2) %>% as.vector()
  cor.param$corr <- sapply(cor.param, function(x) cor.test(df2[[param]],df2[[x]])$estimate)
  cor.param$pvalue <- sapply(cor.param[,1], function(x) cor.test(df2[[param]],df2[[x]])$p.value)
  cor.param <- cor.param[order(cor.param$corr, decreasing = T),]
  
  #plot cor
  df3 <- cor.param %>% filter(pvalue <= 0.1)
  g <- ggplot(data = df3[-1,],aes(x=reorder(predictor,order(df3[-1,]$corr,decreasing=T)),y=corr,label=predictor))+
    geom_col(aes(fill=pvalue))+ 
    geom_text(angle=90,hjust=0,nudge_y = 0.01,size=6)+
    theme_bw(base_size=20)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.grid.minor = theme_line(colour = NA))+
    scale_fill_gradientn(colors = c("skyblue","dodgerblue4","firebrick"),
                         breaks = c(0,0.01,0.05, 0.1),labels = c("0","0.01","0.05","0.1"),limits = c(0,0.1), 
                         guide="legend")+
    ylim(-1.1,1.1)
  if(title == T){
    g <- g + labs(x="",y=name_df,fill="p-value")
    print(g)
  } else if(title == F){
    g <- g + labs(x="",y=param,fill = "p-value")
    print(g)
  }
  
}

# print.cor.signif2 -----
# only prints the significant correlations on the graph
# returns graphs but does not change the list "all.correlations"

print.cor.signif2 <- function(df,param,title = T,y1=-1.1,y2=1.1){
  
  #save name of the dataframe as a character
  df2 <- df
  name_df <- paste(deparse(substitute(df)),param,sep=".")
  
  # correlation table with p values
  cor.param <-data.frame(predictor = rep(NA,length(names(df2))))
  cor.param$predictor <- names(df2) %>% as.vector()
  cor.param$corr <- sapply(cor.param[,1], function(x) cor.test(df2[[param]],df2[[x]])$estimate)
  cor.param$pvalue <- sapply(cor.param[,1], function(x) cor.test(df2[[param]],df2[[x]])$p.value)
  cor.param <- cor.param[order(cor.param$corr, decreasing = T),]
  
  #plot cor
  df3 <- cor.param %>% filter(pvalue <= 0.05)
  g <- ggplot(data = df3[-1,],aes(x=reorder(predictor,order(df3[-1,]$corr,decreasing=T)),y=corr,label=predictor))+
    geom_col(aes(fill=pvalue))+ 
    geom_text(angle=90,hjust=0,nudge_y = 0.01,size=7)+
    theme_bw(base_size=28)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),panel.grid.minor = element_blank())+
    scale_fill_gradientn(colors = c("skyblue","dodgerblue4","firebrick"),
                         breaks = c(0,0.01,0.025, 0.05),labels = c("0","0.01","0.025","0.05"),limits = c(0,0.05), 
                         guide="legend")+
    xlab("")+ylab("param")+
    ylim(y1,y2)
  if(title == T){
    g <- g + labs(x="", y=name_df, fill="p-value")
    print(g)
  } else {
    g <- g + labs(x = "",y = param, title = title)
    print(g)
  }
  
}

# plot.radar ----
# plot radar of the outliers for one parameter -----

plot.radar <- function(outlier_df,param){
  
  library("fmsb")
  library("colormap")
  
  nb_param <- dim(outlier_df)[2]
  nb_sample <- dim(outlier_df)[1]
  nbrow <- ceiling(nb_sample/3)
  outlier_ordered <- outlier_df[order(outlier_df[,param], decreasing = F),]
  outminmax <- rbind(sapply(outlier_ordered, min), sapply(outlier_ordered, max))
  colors_border <- colormap(colormap = colormaps$viridis, nshades = nb_sample, alpha = 1)
  colors_in <- colormap( colormap = colormaps$viridis, nshades = nb_sample, alpha = 0.3)
  par(mfrow = c(nbrow,3), mar = c(0,0,0,0) + 2)
  lapply(1:nb_sample, function (i){radarchart(rbind(outminmax,outlier_ordered[i,]),pfcol = colors_in[i],pcol=colors_border[i], vlcex = 2, cglty = 1, cglcol = "lightgrey" )})

}

# plot.radar2 : radarplot with title -----

plot.radar2 <- function(outlier_df,param){
  outlier_df2 <- outlier_df[,-1]
  nb_param <- dim(outlier_df2)[2]
  nb_sample <- dim(outlier_df)[1]
  nbrow <- ceiling(nb_sample/3)
  outlier_ordered <- outlier_df2[order(outlier_df2[,param], decreasing = F),]
  outminmax <- rbind(sapply(waterchem, min, na.rm = T), sapply(waterchem, max, na.rm = T)) %>% as.data.frame() %>% dplyr::select(names(outlier_df2))
  colors_border <- colormap(colormap = colormaps$viridis, nshades = nb_sample, alpha = 1)
  colors_in <- colormap( colormap = colormaps$viridis, nshades = nb_sample, alpha = 0.3)
  #png(paste(param,"outliers","radarchart","png",sep = "."), width= 800, height = 600)
  par(mfrow = c(nbrow,3), mar = c(0,0,0,0) + 1)
  lapply(1:nb_sample, function (i){radarchart(rbind(outminmax,outlier_ordered[i,]),pfcol = colors_in[i],pcol=colors_border[i], vlcex = 2, cglty = 2, cglcol = "lightgrey", title = outlier_df[i,1] )})
  #dev.off()
}


# myplot -----
# plots comparisons -----
myplot <- function (xvar,yvar,data,lab.x = xvar,lab.y = yvar){
  
  lmxy <- lm(data[,yvar]~data[,xvar]) %>% summary()
  ggplot(data,aes_string(xvar,yvar))+
    geom_point(size=1)+
    geom_text_repel(aes(label = Sample_ID, col = as.character(niva.incub_date)),size = 3,show.legend = F)+
    theme_bw(base_size = 20)+geom_abline(intercept = lmxy$coefficients[1], slope = lmxy$coefficients[2], col = "grey50")+
    scale_color_viridis_d(end = 0.8, direction = -1)+
    labs(x = lab.x, y = lab.y, title = paste("R2 = ",round(lmxy$r.squared,digits = 3),sep = ""))

  }

# fill.na ----
# replaces NA by the mean value of all columns

fill.na <- function(df){
  for (i in names(df)){
    lna <- which(is.na(df[,i]))
    df[lna,i] <- rep(mean(df[,i],na.rm=T),length(lna))
  }
  return(df)
}

# detects na
 not_all_na <- function(x) any(!is.na(x))
 not_any_na <- function(x) all(!is.na(x))
 