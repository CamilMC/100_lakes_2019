library(readxl)
library(dplyr)
library(reshape2)
library(ggplot2)
library(factoextra)
library(viridis)
library(dendextend)
library(circlize)

# load bacteria data

lakes73 <- read.csv("Data/rr100lakes.csv")
seq_asv <- read_xlsx("Data/CBA_100Lakes_bacteria.xlsx")

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
key <- read_xlsx("Data/CBA_100Lakes_bacteria.xlsx",sheet="key") %>% select(c("Lake_ID","Lake_name","DNAsampleID"))
bacteria_lakes_filtered <- dplyr::filter(bacteria_lakes,DNAsampleID != setdiff(bacteria_lakes$DNAsampleID,key$DNAsampleID))
bacteria_73lakes <- merge(key,bacteria_lakes_filtered,by.x="DNAsampleID",by.y="DNAsampleID")
bacteria_73lakes$DNAsampleID <- NULL

# dataframe without ID included
bacteria73_noid <- bacteria_73lakes
rownames(bacteria73_noid) <- bacteria_73lakes$Lake_name
bacteria73_noid$Lake_name <- NULL

# calculate distance of each sample to Langtjern
dist_bacteria <- dist(bacteria73_noid) %>% as.matrix() %>% as.data.frame() %>% select("Langtjern") %>% setNames("dist_bacteria")
dist_bacteria$Sample_ID <- rownames(dist_bacteria)

# makes clusters out of dist
hc <- hclust(dist(bacteria73_noid), method = "ward.D2")

# optimal number of clusters -----
bact73_scaled <- scale(bacteria73_noid) %>% as.data.frame()

# elbow method
fviz_nbclust(bact73_scaled, kmeans, method = "wss") 

#silhouette method
fviz_nbclust(bact73_scaled, kmeans, method = "silhouette")

fviz_nbclust(bact73_scaled, kmeans, method = "gap_stat")

#plot
nb_clusters <- 3
colors_cluster <- viridis(nb_clusters,end = 0.8,direction = -1)

dend_bacteria <- as.dendrogram(hc)
par(mar=c(5,1,1,1))
plot(dend_bacteria)

dend_bact_polar <- dend_bacteria_name %>% color_branches(k=nb_clusters,col=colors_cluster) %>% color_labels(k=nb_clusters,col=colors_cluster)

png("Figures/dendrogram_bact.png",type = "cairo", units = "in", res = 150, width = 8, height = 8)
circlize_dendrogram(dend_bact_polar,dend_trackheight = .8,labels_track_height = 0.3)
dev.off()


# adds cluster to lakes73 dataframe

lakes73$bact_cluster <- cutree(dend_bacteria,nb_clusters)
aov(RR~bact_cluster,data=lakes73) %>% summary()
aov(BdgT~bact_cluster,data=lakes73) %>% summary()
aov(RR/DOC~bact_cluster,data=lakes73) %>% summary()



                                                     