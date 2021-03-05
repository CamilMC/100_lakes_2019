library("scam")
library("DescTools")


# loop

RR <- NULL
pdf("doxdt.pdf")

for (i in 1:length(datalist)) {
  
  d <- datalist[[i]] #select the dataframes of the list one by one
  X <- d[,-which(names(d)==c("Time.Min.","TempC"))] #creates a matrix X containing only the vials (no temperature and no time)
  X <- X[,! apply(X,2,function(x) all(is.na(x)))] #remove columns of vials containing only NA
  X[is.na(X)] <- 0 #removes NA is the other columns
  vials <- names(X) #creates a vector with the names of the selected vials
  t <- d[,"Time.Min."] / 60 #creates a vector with the time values in hours
  m <- vector(length(vials), mode = "list") # creates a matrix with as many columns as vials
  RRi <- data.frame(matrix(NA,nrow=length(vials),ncol=13)) #creates a dataframe to store the results from the fit/rr analysis
  names(RRi) <- c("sample","fit_r2","ox_initial","t_initial","dmax","tmax_h","se","auc","hauc","width","t1_h","t2_h","hwidth") #names of the dataframe
  
  for (j in 1:length(vials)) { #compute the scam fit for all columns in X and stores the fit in m
    m[[j]] <- scam(X[,j] ~ s(t, k=10, bs = "mpd"))
  }
  
  Xs <- sapply(m, predict) #creates matrix with the model
  dXdt <- sapply(m, function(x) { diff(predict(x)) }) # compute derivative of the model
  se.dXdt <- sapply(m, function(x) { derivative.scam(x)$se.d }) #compute standard error of the derivative
  
  #
  
  for (j in 1:length(vials)){
    
    name_sample <- vials[j]  #ID of the sample
    doxdt <- -dXdt[,j]/(60/3600) #column of the model corresponding to the vial (change sign to get the speed of oxygen consumption)
    l <- length(doxdt)
    se.doxdt <- se.dXdt[,j] #corresponding se
    dmax <- max(doxdt) #maximum speed
    xmax <- which.max(doxdt) #abciss of maximum speed
    doxdt1 <- doxdt[1:xmax]
    doxdt2 <- doxdt[xmax:l]
    
    
    x1 <- tail(which(doxdt1 < (dmax/2)),n=1) #abciss of first half-peak point
    if(length(x1)==0){x1 <- xmax} #if x1 is na, the value of xmax is attributed to x1
    t1 <- t[x1] #time in hour corresponding to the abciss
    
    x2 <- xmax -1 + (which(doxdt2 < (dmax/2))[1]) #abciss of second half-peak point
    if(is.na(x2)){x2 <- xmax} #if x2 is na, the value of xmax is attributed to x2
    t2 <- t[x2] #time in hour corresponding to the abciss
    
    A <- AUC(t[-1],doxdt,from = t1, to = t2 ,method = "step") #area under the peak from half-height
    
    if (doxdt[x1] > dmax/2){
      tstart <- t[xmax]
      tstop <- t2
    } else if (doxdt[x2] > dmax/2){
      tstart <- t1
      tstop <- t[xmax]
    } else if (length(x1:xmax) > length(xmax:x2)){
      tstart <- t[xmax]
      tstop <- t2
    } else {
      tstart <- t1
      tstop <- t[xmax]
    }
    
    hA <- AUC(t[-1],doxdt,from = tstart, to = tstop, method = "step")
    
    #plots
    par(mfrow=c(2,1),mar= c(4,4,1,3))
    matplot(x=t, y=X[,j], type="p", pch=".", col="black", xlab="time (h)", ylab="O2 (?M)")
    matplot(x=t, y=Xs[,j], type="l", lty=1, xlab = "", ylab = "", col="red", add=TRUE)
    mtext(paste("R2 =",format(summary(m[[j]])$r.sq,digits=4)), side = 4, line = 1)
    matplot(t[-1], doxdt, type="l", lty=1, col="black", xlab=name_sample, ylab="dO2/dt (?M/h)",ylim=c(min(doxdt-se.doxdt),max(doxdt+se.doxdt)))
    points(t1,doxdt[x1],col="black")
    points(t2,doxdt[x2],col="red")
    abline(v=tstart,col="green")
    abline(v=tstop,col="green")
    abline(h=(dmax/2),col="blue")
    
    
    #fills dataframe with parameters
    RRi$sample[j] <- paste(names(datalist)[i],name_sample,sep="_")
    RRi$fit_r2[j] <- summary(m[[j]])$r.sq
    RRi$ox_initial[j] <- Xs[j][1]
    RRi$t_initial[j] <- t[1]
    RRi$dmax[j] <- dmax
    RRi$tmax_h[j] <- t[xmax]
    RRi$se[j] <- se.doxdt[xmax]
    RRi$auc[j] <- A
    RRi$hauc[j] <- hA
    RRi$width[j] <- t2-t1
    RRi$t1_h[j] <- t1
    RRi$t2_h[j] <- t2
    RRi$hwidth[j] <- tstop-tstart
  }
  
  RR <- rbind(RR,RRi)
  
}

dev.off()

write.csv(RR,"all_RR.csv") #save RR
