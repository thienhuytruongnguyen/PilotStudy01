##Thien 14/06/2021
##Generate comparison plots for Mean, Std Dev, 5th percentile, 95th percentile 
#of monthly total and mean, std dev of no. of wet day
##Detail<<input are observed rainfall data frames and a list contains monthly simulated data frames
##Detail<<require additional functions: pred.interval() , percentile 5th and percentile 95th
monthlytotal_stats_plot<-function
(monthly.obs, ##<< is a dataframe of observed rainfall at a gauged with the following format (column 1: Date; column 2: rainfall amount (mm); column 3: day; column 4:month; column 5: year)
 monthly.pred.list,##<<A list that contain 12 dataframes corresponding for 12 month each dataframe contains the simulated rainfall values for each month with length equal to observed data.frames
 threshold = 0
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
  ylab.plot1 <- c("mean monhtly total rainfall (mm)","std dev monhtly total rainfall (mm)","5th percentile monhtly total rainfall (mm)",
                  "95th percentile monhtly total rainfall (mm)","mean monthly number of wet day","std dev monthly number of wet day")
  for (i in 1:6){
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
    title(ylab=ylab.plot1[i])
    
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


