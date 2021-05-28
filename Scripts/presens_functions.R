#-------------------------------------------------------------------------------
#    IMPORT_OX

# xlsx_name: name of the file to import (file exported from SDR software)
# name_sheet: name of the sheet containing the data (typically "Oxygen Orginal)
#-------------------------------------------------------------------------------

import_ox <- function (xlsx_file, name_sheet){

  ox <- read_excel(path = xlsx_file, sheet = name_sheet)

  if (names(ox[1]) != "Date/Time"){
    colnames(ox) <- unlist(ox[12, ])
    ox <- ox[-(1:12), -(27:30)]
    ox <- ox[ ,-28]

  }
  
  
  names(ox)[grepl("Min",names(ox)) == TRUE] <- "Time.Min."
  names(ox)[grepl("T_internal", names(ox)) == TRUE] <- "TempC"
  ox <- ox[ , -1]
  ox <- as.data.frame(sapply(ox,as.numeric))
  
}

#-------------------------------------------------------------------------------
#    VIALPLOT_LIST

# df_list: a list of the dataframes obtained from the .xlsx containing SDR results (see function "load_data"). MUST HAVE NAMES
# row_col: must be "col" if samples are organised in columns (4 replicates in lines 1 to 6)
#          and "row" if samples are organised in rows (6 replicates in lines A to D)
# y_lim: a vector with the min and max values for the y-axis

#-------------------------------------------------------------------------------

vialplot_list <- function (df_list,row_col,y_lim) {
pdf("all_plots.pdf", width = 11, height = 8)
  
  for (x in 1:length(df_list)) {
    data <- df_list[[x]]

    vials <- names(data)[-which(names(data)==c("Time.Min.","TempC"))]

    if (row_col == "col") {

      for (line_number in c("1","2","3","4","5","6")) {
        
        vial_nb <- grep(line_number,names(data))
        vials <- names(data)[vial_nb]
        
        plottitle <- paste(names(df_list)[x],"_L",as.character(line_number),sep="")
        
        g <- ggplot(data,aes(x=Time.Min./60))+
          geom_line(aes(y=data[,vials[1]],color=as.character(vials[1])))+
          geom_line(aes(y=data[,vials[2]],color=as.character(vials[2])))+
          geom_line(aes(y=data[,vials[3]],color=as.character(vials[3])))+
          geom_line(aes(y=data[,vials[4]],color=as.character(vials[4])))+
          labs(x="time (h)",y=expression(paste("[O2] (",mu,"mol/L)")),title=plottitle)+
          ylim(y_lim)+theme_bw()+theme(legend.title=element_blank())
        print(g)
      }
      
    } else if (row_col == "row"){
      
      for (line_letter in c("A","B","C","D")){
        
        vial_nb <- grep(line_letter,names(data))
        vials <- names(data)[vial_nb]
        
        plottitle <- paste(names(df_list)[x],"_L",as.character(line_letter),sep="")
        
        g <- ggplot(data,aes(x=Time.Min./60))+
          geom_line(aes(y=data[,vials[1]],color=as.character(vials[1])))+
          geom_line(aes(y=data[,vials[2]],color=as.character(vials[2])))+
          geom_line(aes(y=data[,vials[3]],color=as.character(vials[3])))+
          geom_line(aes(y=data[,vials[4]],color=as.character(vials[4])))+
          geom_line(aes(y=data[,vials[5]],color=as.character(vials[5])))+
          geom_line(aes(y=data[,vials[6]],color=as.character(vials[6])))+
          labs(x="time (h)",y=expression(paste("[O2] (",mu,"mol/L)")),title=plottitle)+
          ylim(y_lim)+theme_bw()+theme(legend.title=element_blank())
        print(g)
      }
      
    } else {
      print("Please define row_col argument by 'row' or 'col'")
    }
  }
  dev.off()
}


#-------------------------------------------------------------------------------
#    TEMP_REV_LIST

# datalist: datalist obtained from the .xlsx containing SDR results (see function "import_ox")

#-------------------------------------------------------------------------------
rev.smooth <- function(x){
  rev(cumsum(rev(x)) / 1:length(x))
}

temp_rev_list <- function(datalist){
  list_stab <- list()
  
  for (x in 1:length(datalist)){
    df <- datalist[[x]]
    l <- dim(df)[1] %>% as.numeric()
    temp.revsm <- rev.smooth(df$TempC)
    
    
    if(mean(temp.revsm) > temp.revsm[1]){
      tstab <- which(temp.revsm >= mean(temp.revsm))[1] %>% as.numeric() 
    } else {
      tstab <- which(temp.revsm <= mean(temp.revsm))[1] %>% as.numeric()
    }
    
    dfstab <- df[tstab:l,]
    list_stab[[x]] <- dfstab
    
    # par(mfrow=c(2,1))
    # par(1,1)
    # plot(df$TempC)
    # abline(h=mean(temp.revsm))
    # points(tstab,df$TempC[tstab],col="red")
    # par(1,2)
    # plot(temp.revsm)
    # abline(h=mean(temp.revsm))
    # points(tstab,temp.revsm[tstab],col="red")
  }
  
  names(list_stab) <- names(datalist)
  time_stab <- lapply(list_stab,function(x) {x$Time.Min.[1]/60}) %>% unlist()
  par(mfrow=c(1,1))
  plot(time_stab)
  
  return(list_stab)
}


#-------------------------------------------------------------------------------
#    BLANK_CORRECT_LIST2

# datalist: datalist obtained from the .xlsx containing SDR results (see function "import_ox")
# blank_line : line or column containing the blanks (1 to 6 or "A" to "D")
# standard_vial: Indicate which vial to exclude from the blank line

#-------------------------------------------------------------------------------

blank_correct_list2 <- function(df_list,blank_line,standard_vial) {
  
  df_list_corr <- df_list
  
  for (y in 1:length(df_list)) {
    
    data <- df_list[[y]]
    
    data_corr <- data
    b <- data[, grep(blank_line,names(data))]
    b <- b[,-which(names(b) %in% standard_vial)]
    b <- b[,! apply(b,2,function(x) all(is.na(x)))] #remove columns of vials containing only NA
    
    b_corr <- rowMeans(b)
    b_mean <- mean(b_corr)
    
    for (x in (names(data))){
      if (x != "Time.Min." & x != "TempC"){
        data_corr[x] <- data[x]/(b_corr/b_mean)
      } else {
        data_corr[x] <- data[x]
      }
      
      df_list_corr[[y]] <- data_corr
      
    }
    
    return(df_list_corr)
    
  }
}

#-------------------------------------------------------------------------------
#    BLANK_CORRECT_LIST

# data: dataframe obtained from the .xlsx containing SDR results (see function "load_data")
# blank_line : line or column containing the blanks (1 to 6 or "A" to "D")

#-------------------------------------------------------------------------------

blank_correct_list <- function(df_list,blank_line) {
  
  df_list_corr <- df_list
  
  for (y in 1:length(df_list)) {
    
    data <- df_list[[y]]
    
    data_corr <- data
    b <- data[, grep(blank_line,names(data))]
    b <- b[,! apply(b,2,function(x) all(is.na(x)))] #remove columns of vials containing only NA
    
    b_corr <- rowMeans(b)
    b_mean <- mean(b_corr)
    
    for (x in (names(data))){
      if (x != "Time.Min." & x != "TempC"){
        data_corr[x] <- data[x]/(b_corr/b_mean)
      } else {
        data_corr[x] <- data[x]
      }
      
      df_list_corr[[y]] <- data_corr
      
    }
    
    return(df_list_corr)
    
  }
}

#-------------------------------------------------------------------------------
#    LIST2CSV

# data: dataframe obtained from the .xlsx containing SDR results (see function "load_data")
# blank_line : line or column containing the blanks (1 to 6 or "A" to "D")

#-------------------------------------------------------------------------------

list2csv <- function(list,file_name){
  for (x in 1:length(list)){
    list.x <- list[x]
    write.csv(list.x,file=paste(file_name,x,"csv",sep="."))
  }
}



