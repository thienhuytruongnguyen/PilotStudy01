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
getSimRain <- function(obs.data, rep = 10, mod = "gama"){
rain.data <- format_TimeSeries(obs.data) #formatting the obs data

#Declaring model parameters objects
occur.param <-data.frame()
amount.param <-data.frame()

#Calibrating model parameters
occur.param <- fitMCModel(rain.data)
amount.param <- fitAmountModel(rain.data,mod)

#Simulating
ls.month.sim <- Amount_model(occur.param, amount.param, rep, obs.data)

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
                          warmup
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
getWeatherData <- function(whichStation,start,finish){
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
  
  #Find second shortest distance (closest distance is of point to itself, therefore use second shortest)
  min.d <- apply(d, 1, function(x) order(x, decreasing=F)[2])
  
  #Result
  newdata <- cbind(merg, merg[min.d,], apply(d, 1, function(x) sort(x, decreasing=F)[2]))
  nearestStationNumber<-as.numeric(newdata[1,4])
  nearestStationNumber<-as.character(nearestStationNumber)
  return(nearestStationNumber)
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
