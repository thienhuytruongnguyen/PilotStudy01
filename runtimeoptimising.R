#Optimise run time of SSE_FlowDurationCurve
#
#function component
SSE_OF <- function(theta,
                   indRainDate,
                   paramGR4J,
                   virObsFlow){
  #Passing element in theta to WGEN parameter
  #Occurence model parameters
  startOF <- Sys.time()
  paramMC <- data.frame(matrix(NA,12,2))
  paramMC[,1] <- theta[1:12]; paramMC[,2] <- theta[13:24]
  #Amount model parameters
  paramAmount <- data.frame(matrix(NA,12,2))
  paramAmount[,1] <- theta[25:36]; paramAmount[,2] <- theta[37:48]
  
  #Generate sim rain with given parameters above (a vector)
  startSimRain <- Sys.time()
  simRainRep <- amountModel(occurParam = paramMC,amountParam = paramAmount, indRainDate = indRainDate, rep = 1) #Get Rainfall replicates
  endSimRain <- Sys.time()
  runTimeSimRain <- endSimRain - startSimRain
  print(runTimeSimRain)
  
  #Generate sim flow with sim rain
  #add simRain to paramGR4J options
  paramGR4J[[3]]$Precip <- simRainRep[,1]
  #RunGR4J model with updated sim rain
  startSimFlow <- Sys.time()
  outputGR4J <- runGR4J(paramGR4J)
  endSimFlow <- Sys.time()
  runtimeSimFlow <- endSimFlow - startSimFlow
  print(runtimeSimFlow)
  
  #Get sim flow from output GR4J
  simFlowRep <- outputGR4J$Qsim
  
  #Calculate Exceedance Probability for sim flow and virobs flow
  startSimFDC <- Sys.time()
  simFDC <- getExceedProb(simFlowRep)
  endSimFDC <- Sys.time()
  runTimeSimFDC <- endSimFDC - startSimFDC
  print(runTimeSimFDC)
  
  startVOFDC <- Sys.time()
  virObsFDC <- getExceedProb(virObsFlow)
  endVOFDC <- Sys.time()
  runTimeVOFDC <- endVOFDC - startVOFDC
  print(runTimeVOFDC)
  #Calculate the Sum of square Error
  err <- simFDC$Flow - virObsFDC$Flow
  SSE <- sum(err^2)
  #SSE <- SSE*100
  
  endOF <- Sys.time()
  runTimeOF <- endOF - startOF
  print(runTimeOF)
  return(SSE)
}


iniThetaDataframe <-
  read.csv(file = paste(WD, "/WGENmodel_parameter_eachSite/Site1_410730.csv", sep = ""))
iniTheta <- rep(0,48)
iniTheta[1:12]<-iniThetaDataframe[1:12,1]; iniTheta[13:24]<-iniThetaDataframe[1:12,2]
iniTheta[25:36]<-iniThetaDataframe[1:12,3]; iniTheta[37:48]<-iniThetaDataframe[1:12,4]


sse <- SSE_OF(theta = iniTheta,
                      indRainDate = indRainDate,
                      paramGR4J = paramGR4J,
                      virObsFlow = virObsFlow)




indRainDate <- makeObsDates(RainDat[,1])

amountModel <- function(occurParam,
                        amountParam,
                        rep,
                        indRainDate){
  
  simRainRep <- data.frame(matrix(NA, nrow = indRainDate$nDy, ncol = rep))
  
  for (i in 1:12){
    for (j in 1:rep){
      bin <- MCmodel(length(indRainDate$i.mm[[i]]), occurParam[i,1], occurParam[i,2])
      pred <- rep(0, length(bin))
      randRain <- rgamma(length(bin), amountParam[i,1], amountParam[i,2])
      simRainRep[indRainDate$i.mm[[i]],j] <- bin * randRain
    }
  }
  return(simRainRep)
}

set.seed(68)
simRainRep1 <- amountModel_V2.0(occurParam = paramMC,amountParam = paramAmount, indRainDate = indRainDate, rep = 1)

simRainRep2 <- getSimFlowRep_Opt(theta = theta, paramGR4J = paramGR4J, rep = 100, obsRain = RainDatFormat)
