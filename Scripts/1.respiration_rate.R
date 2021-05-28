###special functions ----
source("Scripts/presens_functions.R")

### 1. Import data ------

library("readxl")
library("dplyr")
library("ggplot2")

bdg_files <- list.files(path="Data/All",pattern="*.xlsx",full.names=TRUE)
names_samples <- sapply(bdg_files, function(x) { substr(x,7,16)})


bdg_list <- lapply(bdg_files,import_ox,name_sheet="Oxygen Orginal")
bdg_30 <- lapply(bdg_list,filter,Time.Min.<=1800) #exclude data over 30 hours

### 2. Remove stabilisation period ------

bdg_stab <- temp_rev_list(bdg_30)


### 3. Blank correction ------

bdg_c <- blank_correct_list2(bdg_stab,1,"D1")
names(bdg_c) <- names_samples

### 4. Plot ------

vialplot_list(bdg_c,"col",c(0,800))

### 5. Calculate rr and returns "all_RR.csv" (! long computing) ------

datalist <- bdg_c
#detach(package="dplyr")
source("Scripts/doxdt.R")

 ### 6.1 Id table ------
library("dplyr")
id_samples <- read.csv("Data/all_samples.csv") # load list of lakes with IDs

id_samples$set <- paste(id_samples$Start_date,id_samples$SDR_number,id_samples$Sample_line,sep = "_")
idtable <- id_samples %>% select(set,Sample_ID,Survey)
write.csv(idtable,"idtable.csv")

### 6.2 Match RR and ID ------

# to run if analysis starts here
RR <- read.csv("Data/all_RR.csv")
idtable <- read.csv("Data/idtable.csv")

# reorganise RR ----
library("stringr")

RR$set <- str_replace(RR$sample,"[A-Z]","")

#check number of samples and missing ones ----
RR$line <- str_sub(RR$set,-1,-1)
noblanks <- filter(RR, line != 1)
nb_noblanks <- length(noblanks$line)/4

# merge RR and ID -> RRt ----
RRt <- merge(RR,idtable,by = "set") # ! does not contains blanks and standards
RRt$tmax <- RRt$tmax_h-RRt$t_initial
RRt <- RRt[,-which(names(RRt) %in% c("X","X.x","X.y"))]
write.csv(RRt,"Data/RRtotal.csv")

### 6.3 Analyse RR ------
library("corrplot")
RRt <- read.csv("Data/RRtotal.csv")
RRt.cor <- cor(RRt[,c("dmax","tmax_h","tmax","auc","hauc","width","hwidth","t1_h","t2_h","ox_initial")])
RRt.cor.plot <- corrplot.mixed(RRt.cor)

### 6.4 Visualization of extracted parameters with PCA -----
library("vctrs")
library("ggplot2")
library("factoextra")

RRpca <- prcomp(RRt[,c("dmax","tmax_h","tmax","auc","hauc","width","hwidth","t1_h","t2_h","ox_initial")],center=TRUE,scale.=TRUE)
#plot(RRpca)
summary(RRpca)
RRpca$rotation %>% as.data.frame() %>% write.csv("4.Respiration_rate/RR.pca.csv")
fviz_pca_var(RRpca,col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel=TRUE,arrowsize=1,labelsize=7)+theme_bw(base_size=20)

fviz_pca_var(RRpca,col.var="contrib", axes = c(1,3),
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel=TRUE,arrowsize=1,labelsize=7)+theme_bw(base_size=20)

# 7. creates a summary table ----

# run if "group_by" does not work ----
# detach(package:plyr)
# library("dplyr")

RRt <- read.csv("Data/RRtotal.csv")
RR1 <- RRt
RR1$Sample_ID <- as.character(RR1$Sample_ID) #convert the sample ID to a character string
RR1$Survey <- as.character(RR1$Survey) #convert the survey name to a character string

RR2 <- RR1 %>% group_by(Sample_ID,Survey) %>%
summarize(dmax=median(dmax),tmax=median(tmax),tmax_h=mean(tmax_h),auc=median(auc),
             hauc=median(hauc),width=median(width),hwidth=median(hwidth),ox_initial = median(ox_initial),
          incub_date = str_sub(set[1],1,6)) #compute mean of rr for each sample

write.csv(RR2,"Data/RRs.csv")

# 8. Creates summary table for RR100
RR100s <- filter(RR2,Survey == "100_lakes")
write.csv(RR100s,"Data/RR100s.csv")

# 9. Creates summary table for RR1000
#RR1000s <- filter(RR2,Survey == "1000_lakes")
#write.csv(RR1000s,"Data/RR1000s.csv")

