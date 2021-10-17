##Format Time Series
format_TimeSeries <- function(x){
library(lubridate)
library(dplyr)
x <- x %>% #Set format for Date column
  mutate(Date_Time = as.Date(Date_Time, format = "%d/%m/%Y")) %>%
  mutate(month = month(Date_Time)) %>%
  mutate(year = year(Date_Time)) %>%
  mutate(day = day(Date_Time))
return (x)
}
##----------------------------------------##
##get Occurrence statistics
occurr_stats <- function(monthlydata){
  
  #Storing only rain data to a data frame for processing
  rain.only.month <- data.frame(matrix(NA,nrow = length(monthlydata)))
  rain.only.month <- monthlydata
  
  P_dw = 0
  P_ww = 0
  Pie_w = 0
  Pie_d = 0
  P_c <- data.frame(matrix(NA,nrow = (length(rain.only.month)-1)))
  #Loop for each site
  
  #print(s)
  n_dw = 0 #Number of dw events at site
  n_ww = 0 #Number of ww events at site
  n_dd = 0 #Number of dd events at site
  n_wd = 0 #Number of wd events at site
  for (i in 2:length(rain.only.month)){
    if (rain.only.month[i] > 0.1 && rain.only.month[i-1] <= 0.1){
      n_dw = n_dw + 1
    }
    if (rain.only.month[i] > 0.1 && rain.only.month[i-1] > 0.1){
      n_ww = n_ww + 1
    }
    if (rain.only.month[i] <= 0.1 && rain.only.month[i-1] <= 0.1){
      n_dd = n_dd + 1
    }
    if (rain.only.month[i] <= 0.1 && rain.only.month[i-1] > 0.1){
      n_wd = n_wd + 1
    }
    
    #Calculate the conditional probability
    P_dw = n_dw/(n_dw+n_dd)
    P_ww = n_ww/(n_ww+n_wd)
    #Unconditional probabilty
    Pie_w = P_dw/(1 + P_dw - P_ww)
    Pie_d = 1 - Pie_w
  }
  sumary <- data.frame(matrix(NA,nrow = 1, ncol = 4))
  sumary[,1] = P_dw ; sumary[,2] = P_ww ; sumary[,3] = Pie_w ; sumary[,4] = Pie_d
  colnames(sumary) <- c("P_DW","P_WW","Pie_W","Pie_D")
  return(sumary)
}
##----------------------------------------##
##fit Occurence Model (2-state Markov Chain model)
fitMCModel<- function
(obs.data ##<<Formatted observed data
){
  ##Declare a dataframe to store the calibrated parameters
  param <- data.frame(matrix(nrow = 12,ncol = 4))
  
  for (i in 1:12){##Loop for each month
    
    x <- occurr_stats(obs.data[obs.data$month==i,2]) #x is a temporary dataframe to store the outcome of function occurr_stats
    param[i,] <- x[1,] #Store the occurrence parameters set for each month
    
  }
  
  return(param)
  
}
##----------------------------------------##
##fit Amount model (gama and exponential distribution)
fitAmountModel<- function
(obs.data,##<<Formated observed data
 mod##<<Choosing which mod. Either: expo or gama
){
  
  if (mod == "expo"){
    
    amount.param <- data.frame(matrix(NA,nrow = 12, ncol = 1))
    
    for (i in 1:12){
      
      x1 <- obs.data[obs.data$month==i,] ##Passing monthly data to a temporary object
      
      x2<-subset(x1, x1$Value!=0) ##Subsetting day with 0 values
      
      ##<<Using optim() to do the MLE of model parameters
      int.theta=1    # Arbitrary initial values
      lower.theta=1e-5    # Lower bounds of parameters
      upper.theta=1e5  # Upper bounds of parameters
      
      opt=optim(fn=expo.loglike, # Function to be optimised
                par=int.theta,      # Initial Parameter Values
                obs=x2[,2],     # Synthetic data passed to obs argument in loglike.multinormal
                control=list(fnscale=-1,trace=2),  # Change from minimization to maximisation problem (fnscale), put some output
                lower=lower.theta,   # Lower bounds of parameters
                upper=upper.theta,    # Upper bounds of parameters
                method="L-BFGS-B"  # Quasi-Newton method that allows Box constraints see ?optim for details
      )
      est.theta=opt[["par"]] # Optimised parameter values
      amount.param[i,] <- est.theta # Store the optimise parameter in the predefined dataframe
    }
    
    return(amount.param)
    
  }
  else if(mod == "gama"){
    
    amount.param <- data.frame(matrix(NA,nrow = 12, ncol = 2))
    
    for (i in 1:12){
      
      x1 <- obs.data[obs.data$month==i,]
      
      x2<-subset(x1, x1$Value!=0)
      
      int.theta=c(0.1,0.1)    # Arbitrary initial values
      lower.theta=c(1e-5,1e-5)    # Lower bounds of parameters
      upper.theta=c(1,1e5)  # Upper bounds of parameters
      
      opt=optim(fn=gamma.loglike, # Function to be optimised
                par=int.theta,      # Initial Parameter Values
                obs=x2[,2],     # Synthetic data passed to obs argument in loglike.multinormal
                control=list(fnscale=-1,trace=2),  # Change from minimization to maximisation problem (fnscale), put some output
                lower=lower.theta,   # Lower bounds of parameters
                upper=upper.theta,    # Upper bounds of parameters
                method="L-BFGS-B"  # Quasi-Newton method that allows Box constraints see ?optim for details
      )
      est.theta=opt[["par"]] # Optimised parameter values
      amount.param[i,] <- est.theta # Store the optimise parameter in the predefined dataframe
      
    }
    
    return(amount.param)
  }
  else{
    print("Incorrect mod selection. Try again with either (expo) or (gama)")
    return()
  }
}
##----------------------------------------##
##Occurence model (WGEN)
monthly_occurr_model <- function(N,PDW,PWW){
  ## Assume P_c = PDW for January
  P_c = PDW
  
  #Assign an dataframe to store the occurrence binary series
  x <- vector(length = N)
  
  #generate a uniform random series U[0,1] to force the occurrence binary series
  U_t <- runif(N,0,1)
  
  for (j in 1:length(U_t)){ #loop for generating the binary occurrence time series
    if (U_t[j] < P_c){ #If statement to sample a 1 or 0 rainfall occurrence from the uniform random series U_t and P_C
      x[j] = 1
    }else{x[j] = 0}
    
    if(x[j] == 1){ #If statement to Update the P_c for the next sampling
      P_c = PWW
    }else{P_c = PDW}
  }
  return(x)
}
##----------------------------------------##
##Amount model (WGEN)
Amount_model <- function
(occur.param,##<<Markov model parameters, PWW and PDW
 amount.param,##<<Amount model parameters, \lambda for exponential distrib or \alpha and \beta for gamma distrib
 rep,##<<Number of replicates
 obs.data##Formatted observed data
){
  rain.data <- obs.data
  ls.monthly.simrain <- list()##List to store sim data for each month
  
  if(length(amount.param[1,]) == 1){ #If exponential
    
    for (i in 1:12){
      print(i)
      ls.monthly.simrain[[i]] <- data.frame(matrix(NA,nrow=length(rain.data[rain.data$month==i,2]),ncol=rep))
      for (j in 1:rep){
        # Create the occurrence binary series
        bin <- monthly_occurr_model(length(rain.data[rain.data$month==i,2]), occur.param[i,1], occur.param[i,2])
        # Declare monthly dataframe to store the SimRain
        pred <- rep(0,length(bin))
        # Loop to generate rainfall amount for rainny day
        for(k in which(bin == 1)){
          
          pred[k] <- rexp(1,amount.param[i,])
        }
        ls.monthly.simrain[[i]][,j] <- pred
      }
    }
    
    return(ls.monthly.simrain)
    
  }
  
  else if (length(amount.param[1,]) == 2){ #If gamma
    for (i in 1:12){
      print(i)
      ls.monthly.simrain[[i]] <- data.frame(matrix(NA,nrow=length(rain.data[rain.data$month==i,2]),ncol=rep))
      for (j in 1:rep){
        # Create the occurrence binary series
        bin <- monthly_occurr_model(length(rain.data[rain.data$month==i,2]), occur.param[i,1], occur.param[i,2])
        # Declare monthly dataframe to store the SimRain
        pred <- rep(0,length(bin))
        # Loop to generate rainfall amount for rainny day
        for(k in which(bin == 1)){
          
          pred[k] <- rgamma(1,shape=amount.param[i,1],rate=1/amount.param[i,2]) # rexp is the exponential distribution sampling function
        }
        ls.monthly.simrain[[i]][,j] <- pred
      }
    }
    return(ls.monthly.simrain)
  }
  else{
    print("check the input parameters data")
    return()
  }
}
##----------------------------------------##
##Log-likelihood functions for gamma and exponential distribition
#Gamma distribution
gamma.loglike <- function(theta,obs.data){
  n = length(obs.data)
  loglike_i = vector(length = n)
  alpha = theta[1]; beta = theta[2]
  for (i in 1: n){
    loglike_i[i] = log(obs.data[i])
  }
  
  loglike_n = (alpha - 1) * sum(loglike_i) - 1/beta * sum(obs.data) - n * alpha *log(beta) - n * log(gamma(alpha))
  #-n * alpha * log(beta) - n * log(gamma(alpha)) - 1/beta * sum(obs.data) + (alpha - 1) * sum(loglike_i)
  if(is.nan(loglike_n)) {
    browser()
  }
  return(loglike_n)
}

#Exponential distribution
expo.loglike <- function(theta,obs.data){
  n = length(obs.data)
  loglike_i = vector(length = n)
  
  for (i in 1:n){
    loglike_i[i] = obs.data[i]
  }
  loglike_n = n*log(theta) - theta*sum(loglike_i)
  if(is.nan(loglike_n)) {
    browser()
  }
  return(loglike_n)
}
##----------------------------------------##
#get SimRain (WGEN model)
getSimRain <- function(obs.data, rep = 100, mod = "gama"){
rain.data <- format_TimeSeries(obs.data) #formatting the obs data

#Declaring model parameters objects
occur.param <-data.frame()
amount.param <-data.frame()

#Calibrating model parameters
occur.param <- fitMCModel(rain.data)
amount.param <- fitAmountModel(rain.data,mod)

#Simulating
ls.month.sim <-Amount_model(occur.param, amount.param, rep, obs.data)

return(ls.month.sim)
}
##----------------------------------------##
##Prediction interval
pred.interval <- function(x,bound=c(0.05,0.95)){
  interval <- vector(length = 3)
  interval[1] <- quantile(x,prob = bound[1])
  interval[2] <- quantile(x,prob = 0.5)
  interval[3] <- quantile(x,prob = bound[2])
  return(interval)
}
##----------------------------------------##
#5th and 95th percentile
percentile5 <- function(x){
  perc <- quantile(x,prob=0.05)
  return(perc)
}

percentile95 <- function(x){
  perc <- quantile(x,prob=0.95)
  return (perc)
}
##----------------------------------------##
##Wet days statistics plot
#Generate comparison plots for Mean, Std Dev, Skew monthly wet days
##Detail<<input are observed rainfall data frames and a list contains monthly simulated data frames
##Detail<<require additional functions: wet_day_stats() and pred.interval()
wetday_monthlystats_plot <- function
(monthly.obs, ##<< is a dataframe of observed rainfall at a gauged with the following format
 ##(column 1: Date; column 2: rainfall amount (mm); column 3: day; column 4:month; column 5: year)
 monthly.pred.list ##<<A list that contain 12 dataframes corresponding for 12 month
 ##each dataframe contains the simulated rainfall values for each month with length equal to observed data.frames
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
  title_plot1 <- c("Monthly wet day amounts (mm): means","Monthly wet day amounts (mm): skew","Monthly wet day amounts (mm): std dev")
  ylab_plot_2 <- c("mean montly wet day amounts (mm)","skew montly wet day amounts (mm)","std dev montly wet day amounts (mm)")
  
  for (i in 1:3){
    ##Type 1
    ##If statements are used to select the x and y limit of the plot
    if (min(pred.stats[[i]])<min(t(obs.monthwetday.stats[i,]))){
      min.range = min(pred.stats[[i]])
    } else {min.range = min(t(obs.monthwetday.stats[i,]))}
    if (max(pred.stats[[i]])>max(t(obs.monthwetday.stats[i,]))){
      max.range = max(pred.stats[[i]])
    } else {max.range = max(t(obs.monthwetday.stats[i,]))}
    ##Box plot function with modifiable whisker range 
    boxplot.ext(pred.stats[[i]],
                ylim=c(min.range,max.range),
                whiskersProb = c(0.05,0.95))
    ##Point the observed data on the plot
    points(t(obs.monthwetday.stats[i,]),col="red",pch=3)
    title(ylab=ylab_plot_2[i])
    
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
##90 Confidence interval
CI90 <- function(x){
  CI <- vector(length=3)
  z = 1.64
  CI[1] = mean(x) - 1.64*(sd(x)/sqrt(length(x)))
  CI[3] = mean(x) + 1.64*(sd(x)/sqrt(length(x)))
  CI[2] = mean(x)
  return(CI)
}
##----------------------------------------##
##Monthly total statistic plot
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
##---------------------------------------##
##calibrate GR4J model
##Format data
getInput <- function(){
  
}
##InputsModel object
InputsModel <- airGR::CreateInputsModel(FUN_MOD = airGR::RunModel_GR4J, DatesR = DATA$DatesR,
                                        Precip = DATA$P, PotEvap = DATA$E)

#RunOptions object
##1.Index Run and WarmUp period
Ind_Run <- seq(which(format(DATA$DatesR, format = "%Y-%m-%d") == "1970-01-01"),
               which(format(DATA$DatesR, format = "%Y-%m-%d") == "2019-02-28"))
Ind_WarmUp <- seq(which(format(DATA$DatesR, format = "%Y-%m-%d") == "1964-01-01"),
                  which(format(DATA$DatesR, format = "%Y-%m-%d") == "1969-12-31"))
##2.Run Option
RunOptions <- airGR::CreateRunOptions(FUN_MOD = airGR::RunModel_GR4J,
                                      InputsModel = InputsModel, IndPeriod_Run = Ind_Run,
                                      IniStates = NULL, IniResLevels = NULL, IndPeriod_WarmUp = Ind_WarmUp)

#InputsCrit object
InputsCrit <- airGR::CreateInputsCrit(FUN_CRIT = airGR::ErrorCrit_NSE, InputsModel = InputsModel, 
                                      RunOptions = RunOptions, VarObs = "Q", Obs = DATA$Q[Ind_Run])

#CalibOptions object
CalibOptions <- airGR::CreateCalibOptions(FUN_MOD = airGR::RunModel_GR4J, FUN_CALIB = airGR::Calibration_Michel)

#CALIBRATION
OutputsCalib <- airGR::Calibration_Michel(InputsModel = InputsModel, RunOptions = RunOptions,
                                          InputsCrit = InputsCrit, CalibOptions = CalibOptions,
                                          FUN_MOD = airGR::RunModel_GR4J)
#Get GR4J parameter
Param <- OutputsCalib$ParamFinalR

#Run GR4J
OutputsModel <- airGR::RunModel_GR4J(InputsModel = InputsModel, RunOptions = RunOptions, Param = Param)

#Plot
plot(OutputsModel, Qobs = DATA$Q[Ind_Run])

#------------------------------Global Optimization------------------------------------------------#
lowerGR4J <- rep(-9.99, times = 4)
upperGR4J <- rep(+9.99, times = 4)

#Differential Evolution
optDE <- DEoptim::DEoptim(fn = OptimGR4J,
                          lower = lowerGR4J, upper = upperGR4J,
                          control = DEoptim::DEoptim.control(NP = 40, trace = 10))
#Particle Swarm
optPSO <- hydroPSO::hydroPSO(fn = OptimGR4J,
                             lower = lowerGR4J, upper = upperGR4J,
                             control = list(write2disk = FALSE, verbose = FALSE))

#-------------------------------Multi-objective optimization---------------------------------------#

InputsCrit_inv <- InputsCrit
InputsCrit_inv$transfo <- "inv"
algo <- "caRamel"
optMO <- caRamel::caRamel(nobj = 2,
                          nvar = 4,
                          minmax = rep(TRUE, 2),
                          bounds = matrix(c(lowerGR4J, upperGR4J), ncol = 2),
                          func = MOptimGR4J,
                          popsize = 100,
                          archsize = 100,
                          maxrun = 15000,
                          prec = rep(1.e-3, 2),
                          carallel = FALSE,
                          graph = FALSE)
param_optMO <- apply(optMO$parameters, MARGIN = 1, FUN = function(x) {
  airGR::TransfoParam_GR4J(x, Direction = "TR")
})
RunOptions$Outputs_Sim <- "Qsim"
run_optMO <- apply(optMO$parameters, MARGIN = 1, FUN = function(x) {
  airGR::RunModel_GR4J(InputsModel = InputsModel,
                       RunOptions = RunOptions,
                       Param = x)
}$Qsim)
run_optMO <- data.frame(run_optMO)

x<-BasinObs
ind_t<-seq(which(format(DATA$DatesR, format = "%Y-%m-%d")=="1984-01-01"),
           which(format(DATA$DatesR, format = "%Y-%m-%d")=="2012-12-31"))
x[,2] <- DATA$P[ind_t]
x[,4] <- DATA$E[ind_t]
x[,6] <- DATA$Q[ind_t]
## loading catchment data
data(L0123001)

## preparation of InputsModel object
InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR4J, DatesR = x$DatesR,
                                 Precip = x$P, PotEvap = x$E)

## calibration period selection
Ind_Run <- seq(which(format(x$DatesR, format = "%Y-%m-%d")=="1990-01-01"),
               which(format(x$DatesR, format = "%Y-%m-%d")=="1999-12-31"))

## preparation of RunOptions object
RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR4J, InputsModel = InputsModel,
                               IndPeriod_Run = Ind_Run)

## calibration criterion: preparation of the InputsCrit object
InputsCrit <- CreateInputsCrit(FUN_CRIT = ErrorCrit_NSE, InputsModel = InputsModel,
                               RunOptions = RunOptions, Obs = x$Qmm[Ind_Run])

## preparation of CalibOptions object
CalibOptions <- CreateCalibOptions(FUN_MOD = RunModel_GR4J, FUN_CALIB = Calibration_Michel)

## calibration
OutputsCalib <- Calibration_Michel(InputsModel = InputsModel, RunOptions = RunOptions,
                                   InputsCrit = InputsCrit, CalibOptions = CalibOptions,
                                   FUN_MOD = RunModel_GR4J)

## simulation
Param <- OutputsCalib$ParamFinalR
OutputsModel <- RunModel_GR4J(InputsModel = InputsModel,
                              RunOptions = RunOptions, Param = Param)

##---------------------------------------##
##Optim GR4J
OptimGR4J <- function(ParamOptim) {
  ## Transformation of the parameter set to real space
  RawParamOptim <- airGR::TransfoParam_GR4J(ParamIn = ParamOptim,
                                            Direction = "TR")
  ## Simulation given a parameter set
  OutputsModel <- airGR::RunModel_GR4J(InputsModel = InputsModel,
                                       RunOptions = RunOptions,
                                       Param = RawParamOptim)
  ## Computation of the value of the performance criteria
  OutputsCrit <- airGR::ErrorCrit_RMSE(InputsCrit = InputsCrit,
                                       OutputsModel = OutputsModel,
                                       verbose = FALSE)
  return(OutputsCrit$CritValue)
}
##Multi-objective optim
MOptimGR4J <- function(i) {
  if (algo == "caRamel") {
    ParamOptim <- x[i, ]
  }
  ## Transformation of the parameter set to real space
  RawParamOptim <- airGR::TransfoParam_GR4J(ParamIn = ParamOptim,
                                            Direction = "TR")
  ## Simulation given a parameter set
  OutputsModel <- airGR::RunModel_GR4J(InputsModel = InputsModel,
                                       RunOptions = RunOptions,
                                       Param = RawParamOptim)
  ## Computation of the value of the performance criteria
  OutputsCrit1 <- airGR::ErrorCrit_KGE(InputsCrit = InputsCrit,
                                       OutputsModel = OutputsModel,
                                       verbose = FALSE)
  ## Computation of the value of the performance criteria
  OutputsCrit2 <- airGR::ErrorCrit_KGE(InputsCrit = InputsCrit_inv,
                                       OutputsModel = OutputsModel,
                                       verbose = FALSE)
  return(c(OutputsCrit1$CritValue, OutputsCrit2$CritValue))
}