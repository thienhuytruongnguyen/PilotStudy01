##Monthly total statistic plot
##Generate comparison plots for Mean, Std Dev, 5th percentile, 95th percentile 
#of monthly total and mean, std dev of no. of wet day
##Detail<<input are observed rainfall data frames and a list contains monthly simulated data frames
##Detail<<require additional functions: pred.interval() , percentile 5th and percentile 95th
monthlytotal_stats_plot<-function
(monthly.obs, ##<< is a dataframe of observed rainfall at a gauged with the following format (column 1: Date; column 2: rainfall amount (mm); column 3: day; column 4:month; column 5: year)
 monthly.pred.list,##<<A list that contain 12 dataframes corresponding for 12 month each dataframe contains the simulated rainfall values for each month with length equal to observed data.frames
 threshold = 0,
 type = "boxplot" #or "errorbar"
){
  ##Passing observed and simulated data from input
  rain.data <- monthly.obs
  ls.monthly.simrain <- monthly.pred.list
  ##Calculate some constants
  rep =  length(ls.monthly.simrain[[1]][1,])##Number of replicate found in the input
  first.year <- min(rain.data$year)##First year of the observed and simulation input
  last.year <- max(rain.data$year)##Last year of the observed and simulation input
  n.year <- last.year - first.year##Total number of years
  ##Declare lists and data frame to be used to store the calculates statistic for each month and each replicates
  month.total <- list()##List contain 12 monthly total dataframes corresponding to 12 month each dataframe has n.year rows and 12 rep columns
  month.wetday <- list()##List contain 12 monthly No. wetday dataframe corresponding to 12 month each dataframe has n.year rows and 12 rep columns
  pred.stats <- list()##List contain 6 dataframes corresponding to 6 stats each dataframe has rep rows and 12 columns
  monthstats.list <- list()##List contain 12 dataframe corresponding to 12 months each data frame has 6 rows and rep columns
  
  #Calculate 6 monthly stats for the observed data------------------------------------------------#
  obs.month.total <- data.frame(matrix(NA,n.year,12))##Monthly total dataframe
  obs.month.total.stats <- data.frame(matrix(NA,6,12))##stats of monthly total dataframe
  month.no.obs.wetday <- data.frame(matrix(NA,n.year,12))##Montly number of wet day
  for (i in 1:12){##Loop for each month
    
    x <- rain.data[rain.data$month==i,]##Passing monthly observed data of all n.year to a temporary object x
    i.year <- first.year ##Year counter start first year end at last year
    for (k in 1:n.year){##Loop for n.year
      x1 <- x[x$year==i.year,]##Passing each monthly of 1 year to a temporary object x
      obs.month.total[k,i]<-sum(x1[,2])##Calculate the monthly total
      #wet.day = 0##Wet day counter for each month of each year
      #for (j in which(x1[,2]>threshold)){#Loop for month (e.g. Feb of 1890 has 28 days, Feb of 1892 has 29 days)
      #couting a wet day if value of that day greater than 0
      wet.day = length(which(x1[,2]>threshold))
      #}
      month.no.obs.wetday[k,i] <- wet.day##Update wet day for the month
      i.year <- i.year + 1##update the year counter and start again for next year
    }
    
  }
  #Calculate the stats for monthly total and no. of wet day
  obs.month.total.stats[1,] <- sapply(obs.month.total,mean)
  obs.month.total.stats[2,] <- sapply(obs.month.total,sd)
  obs.month.total.stats[3,] <- sapply(obs.month.total,percentile5)
  obs.month.total.stats[4,] <- sapply(obs.month.total,percentile95)
  obs.month.total.stats[5,] <- sapply(month.no.obs.wetday,mean)
  obs.month.total.stats[6,] <- sapply(month.no.obs.wetday,sd)
  #------------------------------------------------------------------------------------------------#
  
  #Calculate 6 monthly stats for the simulated data------------------------------------------------#
  ##<<similar to the procedure for the observed input, need 1 more loop for each replicates
  
  for (i in 1:12){##Loop for each month
    
    df <- data.frame(matrix(NA,n.year,rep))##Dataframe for monthly total and no. of wet day n.year values on each rep columns
    month.total[[i]] <- df
    month.wetday[[i]] <- df
    
    df1 <- data.frame(matrix(NA,6,rep))##Dataframe to store the 4 stats of monthly total an 2 stats of no. of wet day
    monthstats.list[[i]] <- df1
    
    month.dat <- ls.monthly.simrain[[i]]##Passing the monthly simulated input to temporary object month.dat
    
    for (j in 1:rep){##Loop for each rep
      x <- rain.data[rain.data$month==i,]#Passing the day format of the observed input to a temporary object x
      x[,2] <- month.dat[,j]##Replace the passing observed value by the simulated value 
      i.year <- first.year#Year counter
      for (k in 1:n.year){#Loop for n.year
        x1 <- x[x$year==i.year,]##Passing each monthly of 1 year to a temporary object x
        month.total[[i]][k,j]<-sum(x1[,2])##Calculate the monthly total
        #wet.day = 0##Wet day counter for each month of each year
        #for (m in which(x1[,2]>threshold)){#Loop for month (e.g. Feb of 1890 has 28 days, Feb of 1892 has 29 days)
        #couting a wet day if value of that day greater than 0
        wet.day = length(which(x1[,2]>threshold))
        #}
        month.wetday[[i]][k,j] <- wet.day##Update wet day for the month
        i.year <- i.year + 1##update the year counter and start again for next year
        
      }
      
    }
    #Calculate the stats for monthly total and no. of wet day
    monthstats.list[[i]][1,] <- sapply(month.total[[i]], mean)
    monthstats.list[[i]][2,] <- sapply(month.total[[i]], sd)
    monthstats.list[[i]][3,] <- sapply(month.total[[i]], percentile5)
    monthstats.list[[i]][4,] <- sapply(month.total[[i]], percentile95)
    monthstats.list[[i]][5,] <- sapply(month.wetday[[i]], mean)
    monthstats.list[[i]][6,] <- sapply(month.wetday[[i]], sd)
  }
  
  ##Gathering and preparing data for plotting  
  for (i in 1:6){##stats loop
    df <- data.frame(matrix(NA,rep,12))##Each stats is a dataframe that has rep rows and 12 columns corresponding for 12 months
    pred.stats[[i]] <- df
    colnames(pred.stats[[i]]) <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
    for (j in 1:12){##Loops for each month
      pred.stats[[i]][,j] <- t(monthstats.list[[j]][i,])##Collect the stats values of each month
    }
  }
  
  ##Start plotting, 2 types of plots: Monthly box plot and Error bar with whisker and range of error bar shown the 90% Prediction interval
  
  title.plot2 <- c("Monthly totals (mm): means","Monthly totals (mm): std dev","Monthly totals (mm): 5th percentile",
                   "Monthly totals (mm): 95th percentile","Monthly no. wet days: means","Monthly no. wet days: std dev")
  title.plot1 <- c("mean monhtly total rainfall (mm)","std dev monhtly total rainfall (mm)","5th percentile monhtly total rainfall (mm)",
                   "95th percentile monhtly total rainfall (mm)","mean monthly number of wet day","std dev monthly number of wet day")
  #par(mfrow=c(2,2))
  for (i in 1:6){
    
    if (type == "boxplot"){
      ##Type 1
      ##If statements are used to select the x and y limit of the plot
      
      if (min(pred.stats[[i]])<min(t(obs.month.total.stats[i,]))){
        min.range = min(pred.stats[[i]])
      } else {min.range = min(t(obs.month.total.stats[i,]))}
      if (max(pred.stats[[i]])>max(t(obs.month.total.stats[i,]))){
        max.range = max(pred.stats[[i]])
      } else {max.range = max(t(obs.month.total.stats[i,]))}
      ##Box plot function with modifiable whisker range 
      boxplot.ext(pred.stats[[i]],
                  ylim=c(min.range,max.range),
                  whiskersProb = c(0.05,0.95))
      ##Point the observed data on the plot
      points(t(obs.month.total.stats[i,]),col="red",pch=3)
      title(main = title.plot1[i])
    } else if (type == "errorbar"){
      
      ##Type 2: Error bar, using arrows() function to mimic error bar
      ## calculate the range of the arrow
      x<- sapply(pred.stats[[i]], pred.interval)
      ## Scatter plot the each 50%percentile point of the sim data to each point of the observed
      plot(x[2,],t(obs.month.total.stats[i,]),pch=3,col="blue",cex=0.5,
           xlim=c(min(x[1,]),max(x[3,])),
           xlab="simulated",ylab="observed")
      ##Plot the arrows with range of 90% prediction interval
      arrows(x0=x[1,],y0=t(obs.month.total.stats[i,]),x1=x[3,],y1=t(obs.month.total.stats[i,]),code=3, angle=90, length=0.0,col = "blue")     
      abline(coef = c(0,1),col="red")
      title(main=title.plot2[i])
    }
  }
}
##---------------------------------------##
##Wet days statistics plot
#Generate comparison plots for Mean, Std Dev, Skew monthly wet days
##Detail<<input are observed rainfall data frames and a list contains monthly simulated data frames
##Detail<<require additional functions: wet_day_stats() and pred.interval()
wetday_monthlystats_plot <- function
(monthly.obs, ##<< is a dataframe of observed rainfall at a gauged with the following format
 ##(column 1: Date; column 2: rainfall amount (mm); column 3: day; column 4:month; column 5: year)
 monthly.pred.list, ##<<A list that contain 12 dataframes corresponding for 12 month
 ##each dataframe contains the simulated rainfall values for each month with length equal to observed data.frames
 type = "boxplot"
){
  
  #require(ggplot2)
  require(moments)
  
  ##Passing observed and simulated data from input
  rain.data <- monthly.obs
  ls.monthly.simrain <- monthly.pred.list
  ##Declare lists and data frame to be used to store the calculates statistic for each month and each replicates
  rep = length(ls.monthly.simrain[[1]][1,]) #Getting number of replicates
  monthstats.list <- list()## List that contains 12 dataframes of monthly statistics for each replicates
  obs.monthwetday.stats <- data.frame(matrix(NA,5,12)) ##This dataframe contains 5 monthly wet day stats for observed data
  colnames(obs.monthwetday.stats)<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")##Names the columns corresponding for each month
  rep.obs <- list()##List contains dataframes that is (rep) replicate of each rows of obs.monthwetday.stats for plotting purpose
  pred.stats <- list() ##List contains dataframes of each month and each replicates stats
  pred.CI.stats <- list()##LIst contains the 90% confidence limit and mean values of simulation sample
  
  ##Calculate simulated monthly wet day statistics 5 statistics: Mean, Skew, Std Dev, Min, Max for each month and for each rep
  for (i in 1:12){## Loop for 12 months
    df <- data.frame(matrix(NA,nrow = 5, ncol = rep)) #Declare property of the statistics dataframe: 
    monthstats.list[[i]] <- df #5 rows for 5 stats, rep columns for number of rep
    
    month.dat <- ls.monthly.simrain[[i]] ##Passing data of each month from the input list
    for(j in 1:rep){ ## Loop for number of replicates
      x<-rain.data[rain.data$month==i,] ## Getting for dates format to work
      x[,2]<-month.dat[,j] ## Replace the obs values with sim values
      monthstats.list[[i]][,j] <- wet_day_stats(x) ## calculating the statistics using the fn wet_day_stats()
    }
  }
  
  ##calculate observed monthly wet day statistics.
  for (i in 1:12){##Loop for 12 months
    obs.monthwetday.stats[,i] <- wet_day_stats(rain.data[rain.data$month==i,]) ##Calculate the statistics using fn wet_day_stats()
    
  }
  
  ##Gathering and preparing data for plotting
  for (i in 1:5){##Loop for each stats
    df1 <- data.frame(matrix(NA,rep,12))##declare property of dataframes (rep) number of rows and 12 columns
    pred.stats[[i]] <- df1 ##rep number of rows and 12 columns
    rep.obs[[i]] <- df1 ##replicates of obs stats with (rep) numbers of identical rows for 12 months
    #pred.CI.stats[[i]] <- data.frame(matrix(NA,3,12))##declare data frames contains upper,lower 90% CI and mean for each month
    for (j in 1:12){##Loop for each month
      pred.stats[[i]][,j] <- t(monthstats.list[[j]][i,])##Passing each stats for each month and each reps to pred.stats dataframe
      #pred.CI.stats[[i]][,j] <- CI90(pred.stats[[i]][,j])##Calculate 90% CI for the pred.stats dataframe using fn CI90()
    }
    colnames(pred.stats[[i]])<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec") ##Names the columns corresponding for each month 
    rep.obs[[i]] <- mapply(rep,obs.monthwetday.stats[i,],rep)##Replicate (rep) number of identical rows for obs stats
  }
  
  ##Start plotting
  title_plot1 <- c("Monthly wet day amounts (mm): means","Monthly wet day amounts: skew","Monthly wet day amounts (mm): std dev")
  title_plot_2 <- c("mean montly wet day amounts (mm)","skew montly wet day amounts","std dev montly wet day amounts (mm)")
  #par(mfrow = c(2,2))
  for (i in 1:3){
    if (type == "boxplot"){
      ##Type 1
      ##If statements are used to select the x and y limit of the plot
      # if (min(pred.stats[[i]])<min(t(obs.monthwetday.stats[i,]))){
      #   min.range = min(pred.stats[[i]])
      # } else {min.range = min(t(obs.monthwetday.stats[i,]))}
      # if (max(pred.stats[[i]])>max(t(obs.monthwetday.stats[i,]))){
      #   max.range = max(pred.stats[[i]])
      # } else {max.range = max(t(obs.monthwetday.stats[i,]))}
      ##Box plot function with modifiable whisker range 
      boxplot.ext(pred.stats[[i]],
                  
                  whiskersProb = c(0.05,0.95))
      ##Point the observed data on the plot
      points(t(obs.monthwetday.stats[i,]),col="red",pch=3)
      title(main = title_plot_2[i])
    } else if (type == "errorbar"){
      ##Type 2: Error bar, using arrows() function to mimic error bar
      ## calculate the range of the arrow
      x<- sapply(pred.stats[[i]], pred.interval)
      ## Scatter plot the each 50%percentile point of the sim data to each point of the observed
      plot(x[2,],t(obs.monthwetday.stats[i,]),pch=3,col="blue",cex=0.5,
           xlim=c(min(x[1,]),max(x[3,])),
           xlab="simulated",ylab="observed")
      ##Plot the arrows with range of 90% prediction interval
      arrows(x0=x[1,],y0=t(obs.monthwetday.stats[i,]),x1=x[3,],y1=t(obs.monthwetday.stats[i,]),code=3, angle=90, length=0.0,col = "blue")     
      abline(coef = c(0,1),col="red")
      title(main=title_plot1[i])
    }
    
  }
}
##----------------------------------------##
##calculate stastistic for wet days
wet_day_stats <- function(x){
  
  require(moments)
  pred.stats <- vector(length = 5)
  
  wet.day.dat <- subset(x,x[,2]!=0)
  
  pred.stats[1] = mean(wet.day.dat[,2])
  pred.stats[2] = skewness(wet.day.dat[,2])
  pred.stats[3] = sd(wet.day.dat[,2])
  pred.stats[4] = min(wet.day.dat[,2])
  pred.stats[5] = max(wet.day.dat[,2])
  
  return(pred.stats)
}
##----------------------------------------##
##Box-plot
boxplot.ext=function(x,whiskersProb=c(0.025,0.975),...){
  # Draws a boxplot with the whiskers at the probability limits, provided by whiskersProb
  # x can be vector or data.frame
  # whiskers Prob values are probabilities for the lower and upper whisker (respectively)
  
  
  x.stats=boxplot(x,plot=FALSE)
  if ( is.data.frame(x) | is.matrix(x)) {
    y=x
    for (j in 1:ncol(x)) {
      y=sort(x[,j])
      if (length(na.omit(y))!=0) {
        x.stats$stats[1,j]=quantile(y,prob=whiskersProb[1])
        x.stats$stats[5,j]=quantile(y,prob=whiskersProb[2])
      }
    }
  }
  else if (is.vector(x)) { # Assume x is a vector
    y=sort(x)
    if (length(na.omit(y))!=0) {
      x.stats$stats[1,1]=quantile(y,prob=whiskersProb[1])
      x.stats$stats[5,1]=quantile(y,prob=whiskersProb[2])
    }
  }
  else
  {
    print("type of x is not supported in boxplot.ext")
    return()
  }
  
  bxp(z=x.stats,...)
  
}
##----------------------------------------##
###Plot virtual observed flow exceedance probability vs sim flow exceedance probability

plotFlowDurationCurve <- function (simFlowRep, 
                                   virObsFlow, 
                                   option = "RepVSVirObs"){ #"RepVSVirObs", "withCILimit"
  
  if (option == "RepVSVirObs"){ ##Not fix for 0 flow yet. will fit if needed in future, to fix see code in "withCILimit" option
    
    #Get exceedance probability for sim flow reps
    simExceedProbRep <- getExceedProbRep(simFlowRep = simFlowRep)
    #Get exceedance probability for virtual observed flow
    virObsExceedProb <- getExceedProb(flow = virObsFlow)
    
    #Plot parameters
    ylab <- "Flow (mm/day)"
    xlab <- "Exceedance Probability"
    #Plot Range
    if (min(simExceedProbRep[[1]]$Flow) < min(virObsExceedProb$Flow)){
      ylower <- min(simExceedProbRep[[1]]$Flow)
    } else {ylower <- min(virObsExceedProb$Flow)}
    
    if (max(simExceedProbRep[[1]]$Flow) > max(virObsExceedProb$Flow)){
      yupper <- max(simExceedProbRep[[1]]$Flow)
    } else {yupper <- max(virObsExceedProb$Flow)}
    
    ylim <- c(ylower,yupper)
    
    #Start plot
    plot(virObsExceedProb, type = "l", log = "y", ylim = ylim, lwd = 2.5, ann = FALSE)
    title(ylab = ylab, xlab = xlab, line = 2.5)
    legend("topright", legend = c("VirObs Flow","Sim Flow", "in Log Scale"),
           col = c("black","blue"), lty = c(1,2,0), lwd = 2, cex = 0.8)
    
    #line replicates
    for (i in 1:ncol(simFlowRep)){
      lines(simExceedProbRep[[i]], lwd = 2,
            col = rgb(red = 0, green = 0, blue = 1, alpha = 0.5),
            lty = 2)
    }
  } else if (option == "withCILimit"){
    
    if (any(virObsFlow==0)==TRUE){
      print("In plotFlowDurationCurve(simFlowRep = simFlowRep, virObsFlow = virObsFlow, option = `RepVSVirObs`) : zeroes detected in 'virObsFlow': some plots in the log space will not be created using all time-steps")
      
      #Get exceedance probability for virtual observed flow
      virObsExceedProb <- getExceedProb(flow = virObsFlow)
      
      #Replace 0 value with NA
      virObsExceedProb$Flow[which(virObsExceedProb$Flow==0)]<-NA
      
      
      if (any(simFlowRep==0)==TRUE){
        print("In plotFlowDurationCurve(simFlowRep = simFlowRep, virObsFlow = virObsFlow, option = `RepVSVirObs`) : zeroes detected in 'simFlowRep': some plots in the log space will not be created using all time-steps")
        
        #Rank simulated flow
        simExceedProbRep <- getExceedProbRep(simFlowRep = simFlowRep)#Rank First time
        
        #Extra 1st ranked to a dataframe
        firstRank <- data.frame(matrix(NA,nrow = length(simExceedProbRep[[1]]$Flow),ncol = ncol(simFlowRep)))
        for (i in 1:ncol(firstRank)){
          firstRank[,i] <- simExceedProbRep[[i]]$Flow
        }
        
        #get probability limit and median 90%
        probLimLower <- apply(firstRank,1,percentile5)
        probLimUpper <- apply(firstRank,1,percentile95)
        probLimMedian <- apply(firstRank,1,percentile50)
        
        #get Exceedance probability for CI limits and median (Rank Second time)
        lowerExceedProb <- getExceedProb(flow = probLimLower)
        upperExceedProb <- getExceedProb(flow = probLimUpper)
        medianExceedProb <- getExceedProb(flow = probLimMedian)
        
        #Replace 0 value with NA
        lowerExceedProb$Flow[which(lowerExceedProb$Flow==0)]<-NA
        upperExceedProb$Flow[which(upperExceedProb$Flow==0)]<-NA
        medianExceedProb$Flow[which(medianExceedProb$Flow==0)]<-NA
        
        newlowerExceedProb <- na.omit(lowerExceedProb)
        newupperExceedProb <- na.omit(upperExceedProb)
        newmedianExceedProb <- na.omit(medianExceedProb)
        
        #Plot parameters
        ylab <- "Flow (mm/day)"
        xlab <- "Exceedance Probability"
        
        #Plot Range
        if (min(lowerExceedProb$Flow, na.rm = TRUE) < min(virObsExceedProb$Flow, na.rm = TRUE)){
          ylower <- min(lowerExceedProb$Flow, na.rm = TRUE)
        }else {ylower <- min(virObsExceedProb$Flow,na.rm = TRUE)}
        
        if (max(upperExceedProb$Flow, na.rm = TRUE) > max(virObsExceedProb$Flow, na.rm = TRUE)){
          yupper <- max(upperExceedProb$Flow, na.rm = TRUE)
        } else {yupper <- max(virObsExceedProb$Flow, na.rm = TRUE)}
        
        ylim <- c(ylower,yupper)
        
      }
      #Start plot
      plot(virObsExceedProb$Exceedance_Probability[!is.na(virObsExceedProb$Flow)],virObsExceedProb$Flow[!is.na(virObsExceedProb$Flow)], log="y", ylim = ylim , type="l", lwd = 2, ann = FALSE)
      title(ylab = ylab, xlab = xlab, line = 2.5)
      legend("topright", legend = c("VirObs Flow","Sim. 90% PL", "Sim. Median", "in Log Scale"),
             col = c("black","red","darkblue"), lty = c(1,2,3,0), lwd = 2, cex = 0.8)
      
      #plot CI region
      xpoly = c(rev(lowerExceedProb$Exceedance_Probability), upperExceedProb$Exceedance_Probability)
      
      ypoly = c(rev(upperExceedProb$Flow), lowerExceedProb$Flow)
      
      polygon(xpoly, ypoly, col = rgb(red = 0.5, green = 0.5, blue = 1, alpha = 0.3), border = NA)
      
      
      #Line CI boundary and median
      lines(newlowerExceedProb, col="red",lwd=1.5, lty=5)                    
      lines(newupperExceedProb, col="red",lwd=1.5, lty=5)
      lines(newmedianExceedProb, col="darkblue",lwd=1.5, lty=3)
    } else{ #Get exceedance probability for virtual observed flow
      virObsExceedProb <- getExceedProb(flow = virObsFlow)
      
      #Rank simulated flow
      simExceedProbRep <- getExceedProbRep(simFlowRep = simFlowRep)#Rank First time
      
      #Extract 1st ranked to a dataframe
      firstRank <- data.frame(matrix(NA,nrow = length(simExceedProbRep[[1]]$Flow),ncol = ncol(simFlowRep)))
      for (i in 1:ncol(firstRank)){
        firstRank[,i] <- simExceedProbRep[[i]]$Flow
      }
      
      
      #get probability limit and median 90%
      probLimLower <- apply(firstRank,1,percentile5)
      probLimUpper <- apply(firstRank,1,percentile95)
      probLimMedian <- apply(firstRank,1,percentile50)
      
      
      #get Exceedance probability for CI limits and median
      lowerExceedProb <- getExceedProb(flow = probLimLower)
      upperExceedProb <- getExceedProb(flow = probLimUpper)
      medianExceedProb <- getExceedProb(flow = probLimMedian)
      
      #Plot parameters
      ylab <- "Flow (mm/day)"
      xlab <- "Exceedance Probability"
      #Plot Range
      if (min(lowerExceedProb$Flow) < min(virObsExceedProb$Flow)){
        ylower <- min(lowerExceedProb$Flow)
      } else {ylower <- min(virObsExceedProb$Flow)}
      
      if (max(upperExceedProb$Flow) > max(virObsExceedProb$Flow)){
        yupper <- max(upperExceedProb$Flow)
      } else {yupper <- max(virObsExceedProb$Flow)}
      
      ylim <- c(ylower,yupper)
      
      #Start plot
      plot(virObsExceedProb, log="y",type="l",ylim = ylim, lwd = 2, ann = FALSE)
      title(ylab = ylab, xlab = xlab, line = 2.5)
      legend("topright", legend = c("VirObs Flow","Sim. 90% PL", "Sim. Median", "in Log Scale"),
             col = c("black","red","darkblue"), lty = c(1,2,3,0), lwd = 1, cex = 0.8)
      #plot CI region
      xpoly = c(rev(lowerExceedProb$Exceedance_Probability), upperExceedProb$Exceedance_Probability)
      ypoly = c(rev(upperExceedProb$Flow), lowerExceedProb$Flow)
      polygon(xpoly, ypoly, col = rgb(red = 0.5, green = 0.5, blue = 1, alpha = 0.3), border = NA)
      
      #Line CI boundary and median
      lines(lowerExceedProb, col="red",lwd=1.5, lty=5)                    
      lines(upperExceedProb, col="red",lwd=1.5, lty=5)
      lines(medianExceedProb, col="darkblue",lwd=1.5, lty=3)
    }
    
    
  }
  
}
#get return interval of annual maxima and plot
compareAnnualMaxima <- function(indObsDate,
                                obs,
                                simRep){
  #Observed Annual maxima
  obsAnnualMaxima <- getAnnualMaxima(indObsDate = indObsDate, value = obs)
  
  #Simulated Annual maxima
  simAnnualMaxima <- data.frame(matrix(NA, nrow = length(obsAnnualMaxima), ncol = ncol(simRep)))
  
  for (i in 1: ncol(simAnnualMaxima)){
    simAnnualMaxima[,i] <- getAnnualMaxima(indObsDate = indObsDate, value = simRep[,i])
  }
  
  #get return interval obsannualMaxima
  annualRetInt_obsAnnualMaxima <- getAnnualRetInt(obsAnnualMaxima)
  
  #get return interval simAnnualMaxima
  annualRetInt_simAnnualMaxima <- getAnnualRetIntRep(simAnnualMaxima)
  
  #get prob limit
  #Extract 1st ranked to a dataframe
  firstRank <- data.frame(matrix(NA,nrow = length(annualRetInt_simAnnualMaxima[[1]]$Depth),ncol = ncol(simRep)))
  for (i in 1:ncol(firstRank)){
    firstRank[,i] <- annualRetInt_simAnnualMaxima[[i]]$Depth
  }
  
  
  #get probability limit and median 90%
  probLimLower <- apply(firstRank,1,percentile5)
  probLimUpper <- apply(firstRank,1,percentile95)
  probLimMedian <- apply(firstRank,1,percentile50)
  
  
  #get Exceedance probability for CI limits and median
  lowerExceedProb <- getAnnualRetInt(dat = probLimLower)
  upperExceedProb <- getAnnualRetInt(dat = probLimUpper)
  medianExceedProb <- getAnnualRetInt(dat = probLimMedian)
  
  #Start plot
  
  #Get range
  if (min(lowerExceedProb) < min(annualRetInt_obsAnnualMaxima)){
    minRange <- min(lowerExceedProb)
  } else {minRange <- min(annualRetInt_obsAnnualMaxima)}
  
  if (max(upperExceedProb) > max(annualRetInt_obsAnnualMaxima)){
    maxRange <- max(upperExceedProb)
  } else {maxRange <- max(annualRetInt_obsAnnualMaxima)}
  
  plot(annualRetInt_obsAnnualMaxima, log="x",pch=4, lwd = 1, ann = FALSE, xaxt ="n", yaxt ="n", ylim = c(minRange, maxRange))
  title(ylab = "Depth (mm)", xlab = "Annual Return Period", line = 2.5)
  
  xticks = c(seq(1,2,0.2),5, 10, 20)
  yticks = seq(0,format(round(max(annualRetInt_obsAnnualMaxima),-1)),40)
  axis(side = 1, at = xticks)
  axis(side = 2, at = yticks)
  abline(h = seq(0, 200, 20), v = xticks, col = "lightgray", lty = 3)
  
  legend("topleft", legend = c("Obs","Sim. 90% PL", "Sim. Median"),
         col = c("black","red","red"), pch = c(4,NA,1), lty = c(0,3,0), lwd = 1, cex = 0.8)
  #Line CI boundary and median
  lines(lowerExceedProb, col="red",lwd=1, lty=3)                    
  lines(upperExceedProb, col="red",lwd=1, lty=3)
  points(medianExceedProb, col="red",cex=0.7, pch = 1)
  
  
}

#Intermittency W/D spell distribution plot
compareDrySpell <- function(monthlyObsRain,
                            monthlySimRain,
                            maxSpell = 10){
  
  #get obs dry Spell
  obsDrySpell <-
    getDrySpell(value = monthlyObsRain, maxSpell = maxSpell)
  
  #get sim dry spell
  simDrySpell <-
    data.frame(matrix(NA, nrow = ncol(monthlySimRain), ncol = maxSpell))
  
  colnames(simDrySpell) <- c(1:10)
  
  for (r in 1:nrow(simDrySpell)){
    simDrySpell[r, ] <-
      getDrySpell(value = monthlySimRain[, r], maxSpell = maxSpell)
  }
  
  #Define ylim
  if (min(simDrySpell) < min(obsDrySpell)){
    minRange <- min(simDrySpell)
  } else {minRange <- min(obsDrySpell)}
  
  if (max(simDrySpell) > max(obsDrySpell)){
    maxRange <- max(simDrySpell)
  } else {maxRange <- max(obsDrySpell)}
  
  #boxplot
  boxplot.ext(
    simDrySpell,
    ylim = c(minRange, maxRange),
    whiskersProb = c(0.05, 0.95)
  )
  title(ylab = "Proportion of Dry Events", xlab = "Dry Spell (Days)")
  
  points(obsDrySpell, pch = 3, col = "red", lwd = 1.5)
  
}

compareWetSpell <- function(monthlyObsRain,
                            monthlySimRain,
                            maxSpell = 10){
  
  #get obs wet Spell
  obsWetSpell <-
    getWetSpell(value = monthlyObsRain, maxSpell = maxSpell)
  
  #get sim wet spell
  simWetSpell <-
    data.frame(matrix(NA, nrow = ncol(monthlySimRain), ncol = maxSpell))
  
  colnames(simWetSpell) <- c(1:10)
  
  for (r in 1:nrow(simWetSpell)){
    simWetSpell[r, ] <-
      getWetSpell(value = monthlySimRain[, r], maxSpell = maxSpell)
  }
  
  #Define ylim
  if (min(simWetSpell) < min(obsWetSpell)){
    minRange <- min(simWetSpell)
  } else {minRange <- min(obsWetSpell)}
  
  if (max(simWetSpell) > max(obsWetSpell)){
    maxRange <- max(simWetSpell)
  } else {maxRange <- max(obsWetSpell)}
  
  #boxplot
  boxplot.ext(
    simWetSpell,
    ylim = c(minRange, maxRange),
    whiskersProb = c(0.05, 0.95)
  )
  title(ylab = "Proportion of Wet Events", xlab = "Wet Spell (Days)")
  
  points(obsWetSpell, pch = 3, col = "red", lwd = 1.5)
  
}