##Format Time Series
format_TimeSeries <- function(x){
require(lubridate)
require(dplyr)
x <- x %>% #Set format for Date column
  mutate(Date_Time = as.Date(Date_Time, format = "%d/%m/%Y")) %>%
  mutate(month = month(Date_Time)) %>%
  mutate(year = year(Date_Time)) %>%
  mutate(day = day(Date_Time))
return (x)
}
##----------------------------------------##
##get Occurrence statistics
occurr_stats <- function(monthlydata,threshold){
  
  #Storing only rain data to a data frame for processing
  #tempMonth <- data.frame(matrix(NA,nrow = length(monthlydata)))
  tempMonth <- monthlydata
  
  P_dw = 0
  P_ww = 0
  Pie_w = 0
  Pie_d = 0
  
  #Loop for each site
  
  #print(s)
  n_dw = 0 #Number of dw events at site
  n_ww = 0 #Number of ww events at site
  n_dd = 0 #Number of dd events at site
  n_wd = 0 #Number of wd events at site
  
  for (i in 2:length(tempMonth)){
    if (tempMonth[i] > threshold && tempMonth[i-1] <= threshold){
      n_dw = n_dw + 1
    }
    if (tempMonth[i] > threshold && tempMonth[i-1] > threshold){
      n_ww = n_ww + 1
    }
    if (tempMonth[i] <= threshold && tempMonth[i-1] <= threshold){
      n_dd = n_dd + 1
    }
    if (tempMonth[i] <= threshold && tempMonth[i-1] > threshold){
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
(obs.data, threshold ##<<Formatted observed data
){
  obs.data <- format_TimeSeries(obs.data)
  ##Declare a dataframe to store the calibrated parameters
  param <- data.frame(matrix(nrow = 12,ncol = 4))
  
  for (i in 1:12){##Loop for each month
    
    tempMonth <- occurr_stats(obs.data[obs.data$month == i, 2], threshold) #x is a temporary dataframe to store the outcome of function occurr_stats
    param[i,] <- tempMonth[1,] #Store the occurrence parameters set for each month
    
  }
  
  return(param)
  
}
##----------------------------------------##
##fit Amount model (gama and exponential distribution)
fitAmountModel<- function
(obs.data,##<<Formated observed data
 mod##<<Choosing which mod. Either: expo or gama
){
  obs.data <- format_TimeSeries(obs.data)
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
      upper.theta=c(1e2,1e5)  # Upper bounds of parameters
      
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

##fit Amount model (gama and exponential distribution)
fitAmountModel_MoM<- function
(obs.data##<<Formated observed data##<<Choosing which mod. Either: expo or gama
){
  obs.data <- format_TimeSeries(obs.data)
  
    
    amount.param <- data.frame(matrix(NA,nrow = 12, ncol = 2))
    
    for (i in 1:12){
      
      monthlyObsData <- obs.data[obs.data$month==i,] ##Passing monthly data to a temporary object
      
      monthlyObsDataWetDay<-subset(monthlyObsData, monthlyObsData$Value!=0) ##Subsetting day with 0 values
      
      ##using method of moment to obtain a and be parameter
      M <- mean(monthlyObsDataWetDay[,2]); V = var(monthlyObsDataWetDay[,2])
      amount.param[i,1] <- (M^2)/V; amount.param[i,2] <- M/V

    }
    return(amount.param)
    }

##----------------------------------------##
##Occurence model (WGEN)
MCmodel <- function(N,PDW,PWW){
  ## Assume P_c = PDW for January
  P_c = PDW
  
  #Assign an dataframe to store the occurrence binary series
  x <- vector(length = N)
  
  #generate a uniform random series U[0,1] to force the occurrence binary series
  #set.seed(68)
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
  #rain.data <- format_TimeSeries(obs.data)
  rain.data <- obs.data
  ls.monthly.simrain <- list()##List to store sim data for each month
  
  if(length(amount.param[1,]) == 1){ #If exponential
    
    for (i in 1:12){
      print(i)
      ls.monthly.simrain[[i]] <- data.frame(matrix(NA,nrow=length(rain.data[rain.data$month==i,2]),ncol=rep))
      for (j in 1:rep){
        # Create the occurrence binary series
        bin <- MCmodel(length(rain.data[rain.data$month==i,2]), occur.param[i,1], occur.param[i,2])
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
        bin <- MCmodel(length(rain.data[rain.data$month==i,2]), occur.param[i,1], occur.param[i,2])
        # Declare monthly dataframe to store the SimRain
        pred <- rep(0,length(bin))
        # locate rain days (indexing)
        indRainDay <- which(bin %in% 1)
        # get number of rain days
        nRainDay <- length(bin[which(bin==1)])
        
        # make rain
        set.seed(68)
        randRain <- rgamma(nRainDay,amount.param[i,1], amount.param[i,2])
        # matching
        
        for (k in 1:nRainDay){
          pred[c(indRainDay)][k] <- randRain[k]
        }
        ls.monthly.simrain[[i]][,j] <- pred
      }
    }
    return(ls.monthly.simrain)
  }
  else{
    print("check the input parameters")
    return()
  }
}
##----------------------------------------##
WGEN_V2.0 <- function(occurParam,
                        amountParam,
                        rep,
                        indRainDate){
  #declare simrain dataframe
  simRainRep <- data.frame(matrix(NA, nrow = indRainDate$nDy, ncol = rep))
  #Loop for each month
  for (i in 1:12){
    #Loop for each replicate
    for (j in 1:rep){
      #Create the occurence binary series
      bin <- MCmodel(length(indRainDate$i.mm[[i]]), occurParam[i,1], occurParam[i,2])
      #make rain ts from gamma distribution
      set.seed(68)
      randRain <- rgamma(length(bin), amountParam[i,1], amountParam[i,2])
      #matching
      simRainRep[indRainDate$i.mm[[i]],j] <- bin * randRain
    }
  }
  return(simRainRep)
}

##----------------------------------------##
WGEN_V3.0 <- function(occurParam,
                             amountParam,
                             rep,
                             indRainDate){
  #declare simrain dataframe
  simRainRep <- data.frame(matrix(NA, nrow = indRainDate$nDy, ncol = rep))
  #Loop for each month
  for (i in 1:12){
    #Loop for each replicate
    for (j in 1:rep){
      #Create the occurence binary series
      bin <- MCmodel(length(indRainDate$i.mm[[i]]), occurParam[i,1], occurParam[i,2])
      #make rain ts from gamma distribution
      #set.seed(68)
      randRain <- rgamma(length(bin[bin==1]), amountParam[i,1], amountParam[i,2])
      
      bin[bin==1] <- randRain
      #matching
      simRainRep[indRainDate$i.mm[[i]],j] <- bin
    }
  }
  return(simRainRep)
}

#---MC model in C++------#
Rcpp::cppFunction('NumericVector MCmodel_C(int n,double PDW, double PWW, NumericVector U_t){
  double P_c = PDW;
  NumericVector day = clone(U_t);
  for (int i = 0; i < n; ++i){
    if (U_t[i] < P_c){
    day[i] = 1;
    } else{
    day[i] = 0;
    }
    if (day[i] == 1){
    P_c = PWW;
    } else{
    P_c = PDW;
    } 
  }
  return day;
}')
##----------------------------------------##
WGEN_V4.0 <- function(occurParam,
                             amountParam,
                             rep,
                             indRainDate){
  #declare simrain dataframe
  simRainRep <- matrix(NA, nrow = indRainDate$nDy, ncol = rep)
  #Loop for each month
  for (i in 1:12){
    #Loop for each replicate
    for (j in 1:rep){
      #Create the occurence binary series
      U_t <- runif(length(indRainDate$i.mm[[i]]),0,1)
      bin <- MCmodel_C(length(U_t), occurParam[i,1], occurParam[i,2], U_t)
      #make rain ts from gamma distribution
      #set.seed(68)
      randRain <- rgamma(length(bin[bin==1]), amountParam[i,1], amountParam[i,2])
      
      bin[bin==1] <- randRain
      #matching
      simRainRep[indRainDate$i.mm[[i]],j] <- bin
    }
  }
  return(simRainRep)
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
getSimRain <- function(obs.data, rep = 10, mod = "gama", option = "MLE", threshold, indRainDate){

  if (option == "MLE"){
    #Declaring model parameters objects
    occur.param <-data.frame()
    amount.param <-data.frame()
    
    #Calibrating model parameters
    occur.param <- fitMCModel(obs.data,threshold)
    amount.param <- fitAmountModel(obs.data,mod)
    
    #Simulating
    simRainRep <- WGEN_V4.0(occur.param, amount.param, rep, indRainDate)
    
    return(list(ls.month.sim, occur.param, amount.param))
  } else if (option == "MoM"){
    #Declaring model parameters objects
    occur.param <-data.frame()
    amount.param <-data.frame()
    
    #Calibrating model parameters
    occur.param <- fitMCModel(obs.data, threshold)
    amount.param <- fitAmountModel_MoM(obs.data)
    
    #Simulating
    simRainRep <- WGEN_V3.0(occur.param, amount.param, rep, indRainDate)
    
    return(list(simRainRep, occur.param, amount.param))
  }

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

percentile50 <- function(x){
  perc <- quantile(x,prob=0.5)
  return(perc)
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
##calibrate GR4J model
##Make input data
makeInputGR4J <- function(
  P,#Rainfall (mm)
  Q,#Flow (ML/day)
  E,#PET (mm)
  Date,
  A #Cacthment area (km^2)
  ){
  #Converting flow unit from ML/day to mm by dividing by catchment area
  Q <- Q/A
  DATA <- data.frame(matrix(NA,nrow = length(RainDat[,1]), ncol = 4))
  colnames(DATA) <- c("DatesR","P","Q","E")
  DATA$DatesR <- Date; DATA$P <- P; DATA$Q <- Q; DATA$E <- E
  
  DATA$DatesR <- strptime(as.character(DATA$DatesR), "%d/%m/%Y")
  DATA$DatesR <- format(DATA$DatesR,"%Y-%m-%d")#format into Y-m-d
  DATA$DatesR <- as.POSIXlt(DATA$DatesR,tz="",format="%Y-%m-%d")#Format date string
  return(DATA)
}
##---------------------------------------##
##calibrate GR4J model
getParamGR4J <- function(inputGR4J,#observed input data
                          start="1970-01-01",#Start of period
                          end="2019-02-28",#End of period
                          warmup,
                         parameter="unknown", knownParam #or "known"
                          ){
  require (airGR)
  require (lubridate)
  
  #get start and end of run and warmup period
  start <- as.Date(start,tryFormats = "%d/%m/%Y"); end <- as.Date(end,tryFormats = "%d/%m/%Y")
  from <- start %m+% months(warmup) ; to <- end
  startWarmUp <- start ; endWarmUp <- from %m-% days(1)
  #set Input model
  InputsModel <- airGR::CreateInputsModel(FUN_MOD = airGR::RunModel_GR4J, DatesR = inputGR4J$DatesR,
                                          Precip = inputGR4J$P, PotEvap = inputGR4J$E)
  #RunOptions object
  ##1.Index Run and WarmUp period
  Ind_Run <- seq(which(inputGR4J$DatesR == from),
                 which(inputGR4J$DatesR == to))
  Ind_WarmUp <- seq(which(inputGR4J$DatesR == startWarmUp),
                    which(inputGR4J$DatesR == endWarmUp))
  ##2.Run Option
  RunOptions <- airGR::CreateRunOptions(FUN_MOD = airGR::RunModel_GR4J,
                                        InputsModel = InputsModel, IndPeriod_Run = Ind_Run,
                                        IniStates = NULL, IniResLevels = NULL, IndPeriod_WarmUp = Ind_WarmUp)
  if (parameter == "unknown"){
    #InputsCrit object
    InputsCrit <- airGR::CreateInputsCrit(FUN_CRIT = airGR::ErrorCrit_NSE, InputsModel = InputsModel, RunOptions = RunOptions, VarObs = "Q", Obs = inputGR4J$Q[Ind_Run])
    
    #CalibOptions object
    CalibOptions <- airGR::CreateCalibOptions(FUN_MOD = airGR::RunModel_GR4J, FUN_CALIB = airGR::Calibration_Michel)
    
    #CALIBRATION
    OutputsCalib <- airGR::Calibration_Michel(InputsModel = InputsModel, RunOptions = RunOptions,
                                              InputsCrit = InputsCrit, CalibOptions = CalibOptions,
                                              FUN_MOD = airGR::RunModel_GR4J)
    #Get GR4J parameter
    Param <- OutputsCalib$ParamFinalR
  } else if (parameter == "known"){
    
    Param <- vector(length = 4); Param[1] <- knownParam[,1]; Param[2] <- knownParam[,2]; Param[3] <- knownParam[,3]; Param[4] <- knownParam[,4]
    
  }

  return(list(Param,Ind_Run,InputsModel,RunOptions))
}
##---------------------------------------##
#run GR4J
runGR4J <- function(paramGR4J){
  require(airGR)
  
  OutputsModel <- airGR::RunModel_GR4J(InputsModel = paramGR4J[[3]], RunOptions = paramGR4J[[4]], Param = paramGR4J[[1]])
  return(OutputsModel)
  # 
  # #Plot
  # plot(OutputsModel, Qobs = inputGR4J$Q[Ind_Run])
}
##---------------------------------------##
## get efficiency criterion
getEffiReport <- function(outputGR4J,inputGR4J,#input data
                          start="1970-01-01",#Start of period
                          end="2019-02-28",#End of period)
                          warmup=12
                          ){
  require(airGR)
  
  #get start and end of run and warmup period
  start <- as.Date(start,tryFormats = "%d/%m/%Y"); end <- as.Date(end,tryFormats = "%d/%m/%Y")
  from <- start %m+% months(warmup) ; to <- end
  startWarmUp <- start ; endWarmUp <- from %m-% days(1)
  
  #set Input model
  InputsModel <- airGR::CreateInputsModel(FUN_MOD = airGR::RunModel_GR4J, DatesR = inputGR4J$DatesR,
                                          Precip = inputGR4J$P, PotEvap = inputGR4J$E)
  #RunOptions object
  ##1.Index Run and WarmUp period
  Ind_Run <- seq(which(inputGR4J$DatesR == from),
                 which(inputGR4J$DatesR == to))
  Ind_WarmUp <- seq(which(inputGR4J$DatesR == startWarmUp),
                    which(inputGR4J$DatesR == endWarmUp))
  
  ##2.Run Option
  RunOptions <- airGR::CreateRunOptions(FUN_MOD = airGR::RunModel_GR4J,
                                        InputsModel = InputsModel, IndPeriod_Run = Ind_Run,
                                        IniStates = NULL, IniResLevels = NULL, IndPeriod_WarmUp = Ind_WarmUp)
  
  #Input criterion
  InputsCrit <- airGR::CreateInputsCrit(FUN_CRIT = airGR::ErrorCrit_NSE, InputsModel = InputsModel, RunOptions = RunOptions, VarObs = "Q", Obs = inputGR4J$Q[Ind_Run])
  
  #Output Criterion report
  OutputsCrit <- ErrorCrit_NSE(InputsCrit = InputsCrit, OutputsModel = outputGR4J)
  
  return(OutputsCrit)
}


##---------------------------------------##
##Getting rain and PET data from URL
getWeatherData <- function(nearestStation,start,finish){
  require(tidyverse)
  require(sf)
  whichStation = list(
    start=start, 
    finish=finish,
    station=nearestStation,
    format="csv",
    comment="RP",
    username="truonghuythien.nguyen@adelaide.edu.au",
    password="silo"
  )
  res <- httr::GET("https://www.longpaddock.qld.gov.au/cgi-bin/silo/PatchedPointDataset.php", query=whichStation)
  silodata <- as.data.frame(read_csv(httr::content(res, as="text")))
  return(silodata)
}
##---------------------------------------##
##Get weather site list
getWeatherSiteList <- function(){
  require(tidyverse)
  require(sf)
  res <- httr::GET("https://www.longpaddock.qld.gov.au/cgi-bin/silo/PatchedPointDataset.php?format=name&nameFrag=_")
  recs <- read_delim(httr::content(res, as="text"), delim = "|"); colnames(recs) <- c('station','name','latitude','longitude','state', 'elevation', 'extra')
   recs <- recs %>%  mutate(latitude = as.numeric(latitude)) %>%
     mutate(longitude = as.numeric(longitude))
  return(recs)
}
##---------------------------------------##
##Get nearest weather station
getNearestWeatherStation <- function(flowSiteList,weatherSiteList,whichSite = 1){
  require(sp)
  require(rgeos)
  #construct spatial classes and perform geo-processing.
  #Read the data and transform them to spatial objects
  flowCoord <- flowSiteList[whichSite,c(1,2,3)]
  siloCoord <- weatherSiteList[,c(1,3,4)]
  merg <- rbind(flowCoord,siloCoord)
  sp.merg <- merg
  coordinates(sp.merg) <- ~longitude+latitude
  
  #calculate pairwise distances between points
  d <- gDistance(sp.merg, byid=T)
  
  #sort nearest silo station
  sort.d <- sort(d[1,2:length(d[,1])], index = TRUE)
  
  #List of 5 nearest SILO station
  nearestStationList <- weatherSiteList[sort.d$ix[1:5],]
  
  #get first silo station info
  infoSilo <- nearestStationList[1,]
  #Result
  nearestStationNumber<-as.numeric(infoSilo[1,1])
  nearestStationNumber<-as.character(nearestStationNumber)
  return(list(nearestStationNumber,infoSilo))
}
##---------------------------------------##
##Get Sim Rain
getSimRainRep <- function(SimRainList){
  rep = length(SimRainList[[1]]) #get number of replicates
  lengthSim <- sum(length(SimRainList[[1]]$X1),length(SimRainList[[2]]$X1),length(SimRainList[[3]]$X1),length(SimRainList[[4]]$X1),length(SimRainList[[5]]$X1),length(SimRainList[[6]]$X1),length(SimRainList[[7]]$X1),length(SimRainList[[8]]$X1),length(SimRainList[[9]]$X1),length(SimRainList[[10]]$X1),length(SimRainList[[11]]$X1),length(SimRainList[[12]]$X1))
  simRain <- data.frame(matrix(NA,nrow = lengthSim, ncol = rep))
  for (i in 1:rep){
    simRain[,i] <- rbind(SimRainList[[1]][i],SimRainList[[2]][i],SimRainList[[3]][i],SimRainList[[4]][i],SimRainList[[5]][i],SimRainList[[6]][i],SimRainList[[7]][i],SimRainList[[8]][i],SimRainList[[9]][i],SimRainList[[10]][i],SimRainList[[11]][i],SimRainList[[12]][i])
  }
  
  return(simRain)
}
##---------------------------------------##
##Flow duration curve
getExceedProb <- function(flow){
  n <- length(flow) #number of events
  P <- data.frame(matrix(NA,nrow = length(flow), ncol = 4))
  colnames(P) <- c("value", "Flow", "Rank", "Exceedance_Probability")
  P[,1] <- flow #get flow value
  P[,2] <- sort(P[,1],decreasing = TRUE) #sort flow value in terms of magnitude
  P[,3] <- 1:nrow(P) #rank
  P[,4] <- P[,3]/(nrow(P)+1) #calculate exceedance probability
  return(P[,c(4,2)])
}
##---------------------------------------##
##Flow duration curve
getExceedProb_V2.0 <- function(flow){

  flow <- flow[order(flow, decreasing = TRUE)]
  Rank <- 1:length(flow)
  #sort flow value in terms of magnitude
  DF <- data.frame(cbind(flow, Rank))
  DF$Prob <- DF$Rank/(length(DF$Rank)+1) #rank
   #calculate exceedance probability
  return(DF[,c(3,1)])
}

##---------------------------------------##
##---------------------------------------##
##Get simflow replicates
getSimFlowRep <- function(simRainRep,paramGR4J){
  require (airGR)
  
  #dataframe to store sim flow rep
  simFlowRep <- data.frame(matrix(NA, nrow = length(paramGR4J[[4]]$IndPeriod_Run), ncol = ncol(simRainRep)))
  
  for (i in 1: ncol(simRainRep)){
    
    paramGR4J[[3]]$Precip <- simRainRep[,i]
  
    outputGR4J <- runGR4J(paramGR4J)
    
    simFlowRep[,i] <- outputGR4J$Qsim
 
  }
  
  return (simFlowRep)
}

##---------------------------------------##
##Get exceedance probability for sim flow replicates

getExceedProbRep <- function(simFlowRep){
  
  repExceedProb <- list()
  
  for (i in 1:ncol(simFlowRep)){
    repExceedProb[[i]] <- getExceedProb(simFlowRep[,i])
  }
  
  return(repExceedProb)
}

##---------------------------------------##
#Manual WGEN
manualWGEN <- function(paramMC,
                       paramAmount,
                       obs.data,
                       rep = 1
){

  SimRainList <- Amount_model(occur.param = paramMC, amount.param = paramAmount, rep = rep, obs.data = obs.data)


  simRainRep <- rbind(SimRainList[[1]],SimRainList[[2]],SimRainList[[3]],SimRainList[[4]],SimRainList[[5]],SimRainList[[6]],SimRainList[[7]],SimRainList[[8]],SimRainList[[9]],SimRainList[[10]],SimRainList[[11]],SimRainList[[12]])


 #return(simRainRep$matrix.NA..nrow...length.rain.data.rain.data.month....i..2....)
  return(simRainRep)
}
##---------------------------------------##
simRaintoFDCModel <- function(paramMC, 
                              paramAmount, obs.rain, 
                              rep = 1,
                              paramGR4J, virObsFlow
){
  #make sim rain
  simRainRep <- manualWGEN(paramMC, paramAmount, obs.rain, rep) #Get Rainfall replicates
  ##Get SimFlow Rep
  simFlowRep <- getSimFlowRep(simRainRep = simRainRep, paramGR4J = paramGR4J)
  #return
  return(simFlowRep)
}
##---------------------------------------##
SSE_FlowDurationCurve_V2.0 <- function(theta,
                              indRainDate,
                              paramGR4J,
                              virObsFlow){
#Passing element in theta to WGEN parameter
  #Occurence model parameters
paramMC <- data.frame(matrix(NA,12,2))
paramMC[,1] <- theta[1:12]; paramMC[,2] <- theta[13:24]
  #Amount model parameters
paramAmount <- data.frame(matrix(NA,12,2))
paramAmount[,1] <- theta[25:36]; paramAmount[,2] <- theta[37:48]

#Generate sim rain with given parameters above (a vector)
  simRainRep <- WGEN_V2.0(occurParam = paramMC,amountParam = paramAmount, indRainDate = indRainDate, rep = 1) #Get Rainfall replicates
#Generate sim flow with sim rain
  #add simRain to paramGR4J options
  paramGR4J[[3]]$Precip <- simRainRep[,1]
  #RunGR4J model with updated sim rain
  outputGR4J <- runGR4J(paramGR4J)
  #Get sim flow from output GR4J
  simFlowRep <- outputGR4J$Qsim
  
#Calculate Exceedance Probability for sim flow and virobs flow
  simFDC <- getExceedProb(simFlowRep)
  virObsFDC <- getExceedProb(virObsFlow)
  
#Calculate the Sum of square Error
  err <- simFDC$Flow - virObsFDC$Flow
  SSE <- sum(err^2)
  #SSE <- SSE*100
  return(SSE)
}
##---------------------------------------##
##
SSE_FlowDurationCurve_V1.0 <- function(theta,
                                       obs.data,
paramGR4J,
virObsFlow){
  #Passing element in theta to WGEN parameter
  #Occurence model parameters
  paramMC <- data.frame(matrix(NA,12,2))
  paramMC[,1] <- theta[1:12]; paramMC[,2] <- theta[13:24]
  #Amount model parameters
  paramAmount <- data.frame(matrix(NA,12,2))
  paramAmount[,1] <- theta[25:36]; paramAmount[,2] <- theta[37:48]
  
  #Generate sim rain with given parameters above (a vector)
  simRainRep <- manualWGEN(paramMC = paramMC, paramAmount = paramAmount, obs.data = obs.data, rep = 1) #Get Rainfall replicates
  #Generate sim flow with sim rain
  #add simRain to paramGR4J options
  paramGR4J[[3]]$Precip <- simRainRep[,1]
  #RunGR4J model with updated sim rain
  outputGR4J <- runGR4J(paramGR4J)
  #Get sim flow from output GR4J
  simFlowRep <- outputGR4J$Qsim
  
  #Calculate Exceedance Probability for sim flow and virobs flow
  simFDC <- getExceedProb(simFlowRep)
  virObsFDC <- getExceedProb(virObsFlow)
  
  #Calculate the Sum of square Error
  err <- simFDC$Flow - virObsFDC$Flow
  SSE <- sum(err^2)
  #SSE <- SSE*100
  return(SSE)
}
##---------------------------------------##
##
SSE_FlowDurationCurve_V3.0 <- function(theta,
indRainDate,
paramGR4J, inputGR4J, runOptionGR4J,
virObsFDC){
  #Passing element in theta to WGEN parameter
  #Occurence model parameters
  paramMC <- data.frame(matrix(NA,12,2))
  paramMC[,1] <- theta[1:12]; paramMC[,2] <- theta[13:24]
  #Amount model parameters
  paramAmount <- data.frame(matrix(NA,12,2))
  paramAmount[,1] <- theta[25:36]; paramAmount[,2] <- theta[37:48]
  
  #Generate sim rain with given parameters above (a vector)
  simRainRep <- WGEN_V3.0(occurParam = paramMC,amountParam = paramAmount, indRainDate = indRainDate, rep = 1) #Get Rainfall replicates
  #Generate sim flow with sim rain
  #add simRain to paramGR4J options
  inputGR4J[[2]] <- simRainRep[,1]
  #RunGR4J model with updated sim rain
  outputGR4J <- airGR::RunModel_GR4J(InputsModel = inputGR4J, RunOptions = runOptionGR4J, Param = paramGR4J)
  #Get sim flow from output GR4J
  simFlowRep <- outputGR4J$Qsim
  
  #Calculate Exceedance Probability for sim flow and virobs flow
  simFDC <- getExceedProb(simFlowRep)
  
  #Calculate the Sum of square Error
  err <- simFDC$Flow - virObsFDC
  SSE <- sum(err^2)
  #SSE <- SSE*100
  return(SSE)
}
##---------------------------------------##
##---------------------------------------##
##
SSE_FlowDurationCurve_V4.0 <- function(theta,
                                       indRainDate,
                                       paramGR4J, inputGR4J, runOptionGR4J,
                                       virObsFDC){
  #Passing element in theta to WGEN parameter
  #Occurence model parameters
  paramMC <- data.frame(matrix(NA,12,2))
  paramMC[,1] <- theta[1:12]; paramMC[,2] <- theta[13:24]
  #Amount model parameters
  paramAmount <- data.frame(matrix(NA,12,2))
  paramAmount[,1] <- theta[25:36]; paramAmount[,2] <- theta[37:48]
  
  #Generate sim rain with given parameters above (a vector)
  simRainRep <- WGEN_V3.0(occurParam = paramMC,amountParam = paramAmount, indRainDate = indRainDate, rep = 1) #Get Rainfall replicates
  #Generate sim flow with sim rain
  #add simRain to paramGR4J options
  inputGR4J[[2]] <- simRainRep[,1]
  #RunGR4J model with updated sim rain
  outputGR4J <- airGR::RunModel_GR4J(InputsModel = inputGR4J, RunOptions = runOptionGR4J, Param = paramGR4J)
  #Get sim flow from output GR4J
  simFlowRep <- outputGR4J$Qsim
  
  #Calculate Exceedance Probability for sim flow and virobs flow
  simFDC <- getExceedProb_V2.0(simFlowRep)
  
  #Calculate the Sum of square Error
  err <- simFDC$flow - virObsFDC
  SSE <- sum(err^2)
  #SSE <- SSE*100
  return(SSE)
}

##---------------------------------------##
##
SSE_FlowDurationCurve_V4.0 <-
  function(theta,
           indRainDate,
           paramGR4J,
           inputGR4J,
           runOptionGR4J,
           virObsFDC) {
    
  
  
  #Passing element in theta to WGEN parameter
  #Occurence model parameters
  paramMC <- matrix(NA,12,2)
  paramMC[,1] <- theta[1:12]; paramMC[,2] <- theta[13:24]
  #Amount model parameters
  paramAmount <- matrix(NA,12,2)
  paramAmount[,1] <- theta[25:36]; paramAmount[,2] <- theta[37:48]
  
  #Generate sim rain with given parameters above (a vector)
  simRainRep <- WGEN_V4.0(occurParam = paramMC,amountParam = paramAmount, indRainDate = indRainDate, rep = 1) #Get Rainfall replicates
  #Generate sim flow with sim rain
  #add simRain to paramGR4J options
  inputGR4J[[2]] <- simRainRep[,1]
  #RunGR4J model with updated sim rain
  outputGR4J <- airGR::RunModel_GR4J(InputsModel = inputGR4J, RunOptions = runOptionGR4J, Param = paramGR4J)
  #Get sim flow from output GR4J
  simFlowRep <- outputGR4J$Qsim
  
  #Calculate Exceedance Probability for sim flow and virobs flow
  simFDC <- getExceedProb_V2.0(simFlowRep)
  
  #Calculate the Sum of square Error
  err <- simFDC$flow - virObsFDC
  SSE <- sum(err^2)
  #SSE <- SSE*100
  return(SSE)
}
##---------------------------------------##
getSimFlowRep_Opt <- function(theta,
                              paramGR4J,
                              rep=100,
                              obsRain){
  
  paramMC <- data.frame(matrix(NA,12,2))
  paramMC[,1] <- theta[1:12]; paramMC[,2] <- theta[13:24]
  paramAmount <- data.frame(matrix(NA,12,2))
  paramAmount[,1] <- theta[25:36]; paramAmount[,2] <- theta[37:48]
  
  indRainDate <- makeObsDates(obsRain[,1])
  simRainRep_Opt <- WGEN_V2.0(paramMC, paramAmount, rep = rep, indRainDate = indRainDate)
  simRainList_Opt <- makeRainList(simRainRep_Opt, indRainDate)
  
  simFlowRep_Opt <-
    getSimFlowRep(simRainRep = simRainRep_Opt, paramGR4J = paramGR4J)
  return(list(simFlowRep_Opt,simRainList_Opt,simRainRep_Opt))
}
# # create index partition of each month, e.g. i.mm[[1]] = c(1,2,3...,31,366,367,...) for Jan; i.mm[[2]]=c(32,33,...)
# # input: dat - observed data vector of dates as string values
# # output is list of length 12 with integer vectors giving indices for relevant months
getMonthlyPartition=function(dat){
  nDy=length(dat) # total number of days in observation record
  dd=rep(0,nDy);mm=rep(0,nDy);yy=rep(0,nDy)  # MAKE EMPTY VECTORS TO STORE DATE INFORMATION SEPARATELY
  for(i in 1:nDy){
    temp=strsplit(as.character(dat[i]),"/")[[1]]  #delimiter is "/" today
    dd[i]=as.integer(temp[1])
    mm[i]=as.integer(temp[2])
    yy[i]=as.integer(temp[3])
  }
  ind=vector(12,mode="list")
  for(m in 1:12) ind[[m]]=which(mm%in%m) #all years
  return(ind)
}

# create index partition of each month, year and month-year
# e.g. i.mm[[1]] = c(1,2,3...,31,366,367,...) for Jan; i.mm[[2]]=c(32,33,...)
#      i.yy[[1]] = c(1,2,3...365)
#      i.ym[[1]][[1]] = c(1,2,...,31) for Jan in year 1
# input: dat - observed data vector of dates as string values
#        CLI - object with climate state information
# output is object with components i.mm i.yy i.ym
makeObsDates=function(dat,CLI=NULL){
  print(paste(Sys.time(),"Start",match.call()[1]))
  nDy=length(dat) # total number of days in observation record
  dd=rep(0,nDy);mm=rep(0,nDy);yy=rep(0,nDy)  # MAKE EMPTY VECTORS TO STORE DATE INFORMATION SEPARATELY
  for(i in 1:nDy){
    temp=strsplit(as.character(dat[i]),"-")[[1]]  #delimiter is "/" today
    dd[i]=as.integer(temp[3])
    mm[i]=as.integer(temp[2])
    yy[i]=as.integer(temp[1])
  }
  # indices for each month
  i.mm=vector(12,mode="list")
  for(m in 1:12) i.mm[[m]]=which(mm%in%m) #all years
  # indices for each year
  stYr=yy[1]
  fnYr=yy[nDy]
  years=stYr:fnYr
  nYr=fnYr-stYr+1
  i.yy=vector(nYr,mode="list")
  for(Y in 1:nYr) i.yy[[Y]]=which(yy%in%years[Y]) # all years
  # indices for each unique year-month
  i.ym=vector(nYr,mode="list") # CREATE YEAR-MONTH INDICES
  for(Y in 1:nYr){
    i.ym[[Y]]=vector(12,mode="list")
    for(m in 1:12){
      i.ym[[Y]][[m]]=i.yy[[Y]][i.yy[[Y]]%in%i.mm[[m]]]
    }
  }
  
  # construction partitions, base "partition" is "all", which means there is no partition at all
  i.part=NULL
  for(Y in 1:nYr) i.part[["all"]]=which(yy%in%years) # all years
  i.mpart=NULL # this is a monthly index lookup of length matching ipart (a lookup of a lookup)
  
  for(m in 1:12)  i.mpart[["all"]][[m]]=which(mm[i.part[["all"]]]%in%m) # all years
  if(!is.null(CLI)){ #there is climate state information to partition the years
    
    for(Y in 1:nYr) i.part[["pos"]]=which(yy%in%CLI$yy.pos) # CLI positive
    for(Y in 1:nYr) i.part[["neg"]]=which(yy%in%CLI$yy.neg) # CLI negative
    
    ###############bug fix here
    for(m in 1:12)  i.mpart[["pos"]][[m]]=i.part[["pos"]][which(mm[i.part[["pos"]]]%in%m)] # CLI positive
    for(m in 1:12)  i.mpart[["neg"]][[m]]=i.part[["neg"]][which(mm[i.part[["neg"]]]%in%m)] # CLI negative
  }
  
  # this is needed for hierarchical model at annual scale, other indices are at daily scale
  # need to rethink the naming convention i.d.ypart (length 129x365) i.m.ypart (length 129x12) i.a.ypart (length 129), similar for i.d.mpart, i.m.mpart
  i.ypart=NULL;i.ypart[["all"]]=1:nYr
  if(!is.null(CLI)){ #there is climate state information to partition the years
    for(Y in 1:nYr) i.ypart[["pos"]]=which(years%in%CLI$yy.pos) # CLI positive - annual timescale
    for(Y in 1:nYr) i.ypart[["neg"]]=which(years%in%CLI$yy.neg) # CLI negative
  }
  
  obsDate=list(i.mm=i.mm,i.yy=i.yy,i.ym=i.ym,nYr=nYr,nDy=nDy,years=years,i.part=i.part,i.mpart=i.mpart,i.ypart=i.ypart,CLI=CLI)
  print(paste(Sys.time(),"Fin",match.call()[1]))
  return(obsDate)
}

# Create indices based on dates for observed data
getDates=function(dat){
  nDy=length(dat) # total number of days in observation record
  dd=rep(0,nDy);mm=rep(0,nDy);yy=rep(0,nDy)  # MAKE EMPTY VECTORS TO STORE DATE INFORMATION SEPARATELY
  for(i in 1:nDy){
    temp=strsplit(as.character(dat[i]),"/")[[1]]  #delimiter is "/" today
    dd[i]=as.integer(temp[1])
    mm[i]=as.integer(temp[2])
    yy[i]=as.integer(temp[3])
  }
  jul=rep(0,nDy);julf=rep(0,nDy) # julian day and julian day fraction
  k=0
  if(isLeap(yy[1])){plusLeapDay=1}else{plusLeapDay=0}
  jul[1]=k;julf[1]=k
  for(i in 2:nDy){
    if(yy[i]!=yy[i-1]){# new year so reset counter
      k=0
      if(isLeap(yy[i])){plusLeapDay=1}else{plusLeapDay=0}
    }else{
      k=k+1
    }
    jul[i]=k
    julf[i]=k/(365+plusLeapDay)
  }
  return(list(dd=dd,mm=mm,yy=yy,jul=jul,julf=julf,nDy=nDy))
}
##-------------------------------------------------------------------#
##Get annual return period
getAnnualRetInt <- function(dat){
  n <- length(dat) #number of events
  P <- data.frame(matrix(NA,nrow = length(dat), ncol = 4))
  colnames(P) <- c("value", "Depth", "Rank", "Annual_Return_Period")
  P[,1] <- dat#get dat value
  P[,2] <- sort(P[,1],decreasing = TRUE) #sort dat value in terms of magnitude
  P[,3] <- 1:nrow(P) #rank
  P[,4] <- (nrow(P)+1)/P[,3] #calculate return interval
  return(P[,c(4,2)])
}

##GEt annual return period rep
getAnnualRetIntRep <- function(simDatRep){
  
  repRetInt <- list()
  
  for (i in 1:ncol(simDatRep)){
    repRetInt[[i]] <- getAnnualRetInt(simDatRep[,i])
  }
  
  return(repRetInt)
}

#---------------------------------#
##Get annual maxima and plot
getAnnualMaxima <- function(indObsDate,
                            value){
  
  annualMaxima <- rep(0,length(indObsDate$i.yy))
  for (i in 1:length(annualMaxima)){
    annualMaxima[i] <- max(value[indObsDate$i.yy[[i]]])
  }
  return(annualMaxima)
}



getSSE <- function(obs, sim){
  err <- data.frame(matrix(NA, nrow = nrow(sim), ncol = ncol(sim)))
  SSE <- rep(0,ncol(sim))
  for (i in 1:ncol(err)){
    err[,i] <- sim[,i] - obs
    SSE[i] <- sum(err[,i]^2)
  }
  SSE <- mean(SSE)
  return(SSE)
} 

##RMSE----#
getRMSE <- function(obs, sim){
  err <- data.frame(matrix(NA, nrow = nrow(sim), ncol = ncol(sim)))
  RMSE <- rep(0,ncol(sim))
  for (i in 1:ncol(err)){
    err[,i] <- sim[,i] - obs
    RMSE[i] <- (sum(err[,i]^2)/length(obs))^(1/2)
  }
  RMSE <- mean(RMSE)
  return(RMSE)
}

##NSE----#
getNSE <- function(obs, sim){
  
  squareErr <- data.frame(matrix(NA, nrow = nrow(sim), ncol = ncol(sim))) #A matrix with rows = number variables, columns = number of replicates
  squareMeanErr <- (obs-mean(obs))^2 #A vector of length = number of variables
  NSE <- rep(0,ncol(sim)) # A vector of length = number of replicates
  
  for (i in 1:ncol(sim)){ #for each replicates calculate the NSE 
    squareErr[,i] <- (sim[,i] - obs)^2
    NSE[i] <- 1 - (sum(squareErr[,i])/sum(squareMeanErr))
  }
  
  NSE <- mean(NSE)#average the NSE of all replicates
  return(NSE)
}

##get NSE for flow duration curve----#
getNSE_FDC <- function(obs,sim){
  
  #Get exceedance probability
  obsExceedProb <- getExceedProb(obs)
  simExceedProb <- getExceedProbRep(sim)
  
  simExceedProbDF <- data.frame(matrix(NA,nrow = length(simExceedProb[[1]]$Flow),ncol = ncol(sim)))
  for (i in 1:ncol(simExceedProbDF)){
    simExceedProbDF[,i] <- simExceedProb[[i]]$Flow
  }
  
  #Calculate NSE
  NSE <- getNSE(obs = obsExceedProb$Flow, sim = simExceedProbDF)
}

##-----------------------Spell Distribution--------------------------------#
getDrySpell <- function(value,
                        maxSpell = 10){
  
  #make binary series of wet and dry day
  bin <- rep(0, length(value))
  
  for (i in which(value > 0)){
    bin[i] <- 1
  }
  
  #calculate wet and dry spells
  spell <- rle(bin)
  
  #get dry spell
  drySpell <- spell$lengths[which(spell$values==0)]
  
  #Totatl dry events
  totalDryEvent <- length(drySpell)
  
  #Calculate dry spell distribution
  drySpellDist <- rep(0, maxSpell)
  
  for (i in 1:maxSpell){
    drySpellDist[i] <- length(drySpell[drySpell==i])/totalDryEvent
  }
  
  return(drySpellDist)
}

##-------------------------------------------------------------------#
getWetSpell <- function(value,
                        maxSpell = 10){
  
  #make binary series of wet and dry day
  bin <- rep(0, length(value))
  
  for (i in which(value > 0)){
    bin[i] <- 1
  }
  
  #calculate wet and dry spells
  spell <- rle(bin)
  
  #get dry spell
  wetSpell <- spell$lengths[which(spell$values==1)]
  
  #Totatl dry events
  totalWetEvent <- length(wetSpell)
  
  #Calculate dry spell distribution
  wetSpellDist <- rep(0, maxSpell)
  
  for (i in 1:maxSpell){
    wetSpellDist[i] <- length(wetSpell[wetSpell==i])/totalWetEvent
  }
  
  return(wetSpellDist)

}  

#Random shuffle a time series
randomShuffle <- function(x,n){
  for (i in 1:n ) {
    x = x[.Internal(sample(length(x),length(x), FALSE, NULL))]
  }
  return(x)
} 
##-------------------------------------------------------------------#

isLeap.vec=function(y){
  y%%4==0& (!(y%%100==0))|(y%%400==0)
}
##-------------------------------------------------------------------#

makeSimDates=function(nYr=10000){
  print(paste(Sys.time(),"Start",match.call()[1]))
  if(nYr>100) print("May take a few minutes ... time for a coffee")
  if(nYr>1000) print("May take a few minutes ... actually, time for lunchbreak")
  print("Creating dates")
  dttmp=seq(as.Date("01/01/2000",format="%d/%m/%Y"), as.Date(paste0("31/12/",2050-1),format="%d/%m/%Y"), "days")
  nsDy=length(dttmp)
  print("Getting date components")
  sdd=as.POSIXlt(dttmp)$mday
  sjul=as.POSIXlt(dttmp)$yday # julian day
  smm=as.POSIXlt(dttmp)$mon+1
  syy=as.POSIXlt(dttmp)$year+1900
  
  sjulf=sjul/(365+as.integer(isLeap.vec(syy))) # julian fraction
  
  print("Creating monthly index")
  j.mm=NULL # CREATE MONTHLY INDICES
  for(m in 1:12) j.mm[[m]]=which(smm%in%m) #all years
  
  # indices for each year
  print("Creating yearly index")
  stYr=syy[1]
  fnYr=syy[nsDy]
  years=stYr:fnYr
  
  j.yy=vector(nYr,mode="list")
  for(Y in 1:nYr) j.yy[[Y]]=which(syy%in%years[Y]) # all years
  
  # indices for each unique year-month
  print("Creating year-month index")
  j.ym=vector(nYr,mode="list") # CREATE YEAR-MONTH INDICES
  for(Y in 1:nYr){
    j.ym[[Y]]=vector(12,mode="list")
    for(m in 1:12){
      j.ym[[Y]][[m]]=j.yy[[Y]][j.yy[[Y]]%in%j.mm[[m]]]
    }
  }
  
  # Now make a chunk at length 400 years - the leap year pattern repeats every 400 years so is convenient for simulation
  # cannot use earlier values since will not work when nYr<400
  nYr.chunk=400 # to simulate 400 years at a time - needed to avoid using up too much memory during long simulation
  # do not vary this parameter because although 100 years of 1000 year chunks might give better performance, extra specialised code will be needed to correct for leap years
  print("Creating quad-century leap year index for simulation chunking")
  dttmp400=seq(as.Date("01/01/0000",format="%d/%m/%Y"), as.Date(paste0("31/12/",nYr.chunk-1),format="%d/%m/%Y"), "days")
  nsDy400=length(dttmp400)
  sdd400=as.POSIXlt(dttmp400)$mday
  sjul400=as.POSIXlt(dttmp400)$yday # julian day
  smm400=as.POSIXlt(dttmp400)$mon+1
  syy400=as.POSIXlt(dttmp400)$year+1900
  
  sjulf400=sjul400/(365+as.integer(isLeap.vec(syy400))) # julian fraction
  
  j.mm400=NULL # CREATE MONTHLY INDICES
  for(m in 1:12) j.mm400[[m]]=which(smm400%in%m) #all years
  
  print("Creating ascii format")
  predatefmt=format(dttmp,format="%d/%m/%Y") # create format for output strings, needed to write ascii files
  predatefmt400=format(dttmp400,format="%d/%m/%Y") # create format for output strings, needed to write ascii files
  
  print(paste(Sys.time(),"Fin",match.call()[1]))
  return(list(j.mm=j.mm,j.yy=j.yy,j.ym=j.ym,years=years,nYr=nYr,str=predatefmt,nsDy=nsDy,smm=smm,sdd=sdd,syy=syy,sjul=sjul,sjulf=sjulf,
              j.mm400=j.mm400,str400=predatefmt400,nsDy400=nsDy400,smm400=smm400,sdd400=sdd400,syy400=syy400,sjul400=sjul400,sjulf400=sjulf400,
              nYr.chunk=nYr.chunk))

}

##---------------------#
getMean3dayTotal <- function(value,#Rainfall timeseries
                             indObsDate){#day index of timeseries
  mean3dayTotal = rep(0,12)
  for (i in 1:12){#Loop for each month
    threeDayTotal <- zoo::rollsum(value[indObsDate$i.mm[[i]]],3) #Get sum 3 day for each month
    mean3dayTotal[i] <- mean(threeDayTotal) #get average 3 day total for each month
  }
  return(mean3dayTotal)
}

##------------------#
getMean5dayTotal <- function(value,
                             indObsDate){
  mean5dayTotal = rep(0,12)
  for (i in 1:12){
    fiveDayTotal <- zoo::rollsum(value[indObsDate$i.mm[[i]]],5)
    mean5dayTotal[i] <- mean(fiveDayTotal)
  }
  return(mean5dayTotal)
}
##---------------------#
makeRainList <- function(simRainRep,
                         indRainDate){
  simRainList <- list()
  
  for (i in 1:12){
    
    tempMatrix <-
      data.frame(matrix(
        NA,
        nrow = nrow(simRainRep[indRainDate$i.mm[[i]],]),
        ncol = ncol(simRainRep)
      ))
    
    for (j in 1:ncol(simRainRep)){
      tempMatrix[,j] <- simRainRep[indRainDate$i.mm[[i]],j]
    }
    simRainList[[i]] <- tempMatrix
  }
  return (simRainList)
}