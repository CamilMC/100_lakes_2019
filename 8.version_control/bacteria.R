# load bacteria data
seq_asv <- read_xlsx("5.100_lakes/bacteria_100Lakes.xlsx")

# vector with class names
bacteria_class <- unique(seq_asv$Class)

# selects only the columns for water samples
seq_asv_class <- seq_asv %>% select(grep("L",names(seq_asv)))

#adds class names
seq_asv_class$Class <- seq_asv$Class
 
#adds the counts by class 
seq_asv_grouped <- seq_asv_class %>% group_by(Class) %>% summarise_all(list(sum))

# transposes the matrix
seq_asv_t <- t(seq_asv_grouped) %>% as.data.frame() %>% setNames(seq_asv_grouped$Class)
seq_asv_t2 <- seq_asv_t[-1,] %>% sapply(as.numeric) %>% as.data.frame()
rownames(seq_asv_t2) <- rownames(seq_asv_t[-1,])
# checking there are no empty columns
which(colSums(seq_asv_t2)==0)

# prepares for match with Lake_ID
bacteria_lakes <- seq_asv_t2
bacteria_lakes$DNAsampleID <- NA
bacteria_lakes$DNAsampleID <- rownames(seq_asv_t2)

#matches with Lake_ID
key <- read_xlsx("5.100_lakes/bacteria_100Lakes.xlsx",sheet="key") %>% select(c("Lake_ID","DNAsampleID"))
bacteria_lakes_filtered <- dplyr::filter(bacteria_lakes,DNAsampleID != setdiff(bacteria_lakes$DNAsampleID,key$DNAsampleID))
bacteria_73lakes <- merge(key,bacteria_lakes_filtered,by.x="DNAsampleID",by.y="DNAsampleID")
bacteria_73lakes$DNAsampleID <- NULL

# Long format
bacteria73 <- melt(bacteria_73lakes,id.vars="Lake_ID",variable.name = "Class") %>% filter(value > 10)

# plot proportional bars
ggplot(bacteria73,aes(y=as.character(Lake_ID),x=value,fill=Class))+geom_col(position = "fill", show.legend = F)+
  theme(legend.position = NULL)+scale_fill_viridis_d(end=0.8,direction = -1)+theme_light()

# plot bars with count
ggplot(bacteria73,aes(y=as.character(Lake_ID),x=value,fill=Class))+geom_col(show.legend = F)+
  theme(legend.position = NULL)+scale_fill_viridis_d(end=0.8,direction = -1)+theme_light()


# dataframe without ID included
bacteria73_noid <- bacteria_73lakes
rownames(bacteria73_noid) <- bacteria_73lakes$Lake_ID
bacteria73_noid$Lake_ID <- NULL

# calculate distance of each sample to Langtjern
dist_bacteria <- dist(bacteria73_noid) %>% as.matrix() %>% as.data.frame() %>% select("13763") %>% setNames("dist_bacteria")
dist_bacteria$Sample_ID <- rownames(dist_bacteria)

# adds the distance to lakes70 dataframe
lakes70_bact <- merge(lakes70,dist_bacteria,by.x="Sample_ID",by.y = "Sample_ID")
lakes70$dist_bact <- lakes70_bact$dist_bacteria

# test difference of distance between clusters 
ggplot(lakes70_bact)+geom_boxplot(aes(x=cluster,y=dist_bacteria,group=cluster))

# reorder columns depending on abundance 
class_abundance <- colSums(bacteria73_noid) %>% as.data.frame()
class_abundance$order <- order(class_abundance,decreasing = T)
bacteria73_noid <- bacteria73_noid %>% select(class_abundance_order)

# explore class abundance
mean_ab <- mean(class_abundance)
sd_ab <- sd(class_abundance)
which(class_abundance < (mean_ab-sd_ab))
quantile(class_abundance,probs = seq(0,1,0.10))
boxplot(class_abundance)

# bacteria73_10classes_long
bacteria73_sorted <- bacteria73_noid
bacteria73_sorted$Lake_ID <- rownames(bacteria73_noid)
bacteria73_sorted <- bacteria73_sorted %>% select("Lake_ID",everything())
bacteria73_sorted_long <- melt(bacteria73_sorted[,1:11],id.vars="Lake_ID",variable.name = "Class")

# plot proportional bars
ggplot(bacteria73_sorted_long,aes(y=as.character(Lake_ID),x=value,fill=Class))+geom_col(position = "fill", show.legend = T)+
  scale_fill_viridis_d()+theme_light()+labs(x="Proportion of the 10 most represented bacterial classes",y="Lake ID")
# plot bars with count
ggplot(bacteria73_sorted_long,aes(y=as.character(Lake_ID),x=value,fill=Class))+geom_col(show.legend = T)+
  scale_fill_viridis_d(end=0.8,direction = -1)+theme_light()+labs(x="Count of the 10 most represented bacterial classes",y="Lake ID")

#checking abundance per sample
ggplot(filter(bacteria73_sorted_long,Lake_ID == 13763))+geom_col(aes(x=Class,y=value))


bact13763 <- bacteria73 %>% filter(Lake_ID == 13763) %>% arrange(desc(value))
ggplot(bact13763,aes(x=Class,y=value,fill=value))+geom_col(show.legend = F)+theme_light()+theme(axis.text.x = element_text(angle=90))+
  scale_x_discrete(limits=bact13763$Class)+scale_fill_viridis_c(direction = -1)


x=reorder(predictor,order(df3[-1,]$corr,decreasing=T))

                                                     