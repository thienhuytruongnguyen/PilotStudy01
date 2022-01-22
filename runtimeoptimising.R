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
############################################################
set.seed(68)
simRainRep1 <- amountModel_V2.0(occurParam = paramMC,amountParam = paramAmount, indRainDate = indRainDate, rep = 1)

simRainRep2 <- getSimFlowRep_Opt(theta = theta, paramGR4J = paramGR4J, rep = 100, obsRain = RainDatFormat)


simRainRep3 <- amountModel_V3.0(occurParam = paramMC,amountParam = paramAmount, indRainDate = indRainDate, rep = 1)

############################################################

SSE <- function(theta,
                indRainDate,
                paramGR4J,
                virObsFlow){
  #s = Sys.time()
  #Passing element in theta to WGEN parameter
  #Occurence model parameters
  #s1 = Sys.time()
  paramMC <- data.frame(matrix(NA,12,2))
  paramMC[,1] <- theta[1:12]; paramMC[,2] <- theta[13:24]
  #e1 = Sys.time()
  #print(e1-s1)
  #Amount model parameters
  #s2 = Sys.time()
  paramAmount <- data.frame(matrix(NA,12,2))
  paramAmount[,1] <- theta[25:36]; paramAmount[,2] <- theta[37:48]
  #e2 = Sys.time()
  #print(e2-s2)
  #Generate sim rain with given parameters above (a vector)
  #s3 = Sys.time()
  simRainRep <- amountModel_V3.0(occurParam = paramMC,amountParam = paramAmount, indRainDate = indRainDate, rep = 1)
  
  #e3 = Sys.time()
  #print(e3-s3)
  #Get Rainfall replicates
  #Generate sim flow with sim rain
  #add simRain to paramGR4J options
  #s4 = Sys.time()
  paramGR4J[[3]]$Precip <- simRainRep[,1]
  #e4 = Sys.time()
  #print(e4 - s4)
  #RunGR4J model with updated sim rain
  #s5 = Sys.time()
  outputGR4J <- runGR4J(paramGR4J)
  #e5 = Sys.time()
  #print(e5-s5)
  #Get sim flow from output GR4J
  #s6 = Sys.time()
  simFlowRep <- outputGR4J$Qsim
  #e6 = Sys.time()
  #print(e6-s6)
  #Calculate Exceedance Probability for sim flow and virobs flow
  #s7 = Sys.time()
  simFDC <- getExceedProb(simFlowRep)
  #e7 = Sys.time()
  #print(e7-s7)
  
  #s8 = Sys.time()
  virObsFDC <- getExceedProb(virObsFlow)
  #e8 = Sys.time()
  #print(e8-s8)
  #Calculate the Sum of square Error
  #s9 = Sys.time()
  err <- simFDC$Flow - virObsFDC$Flow
  SSE <- sum(err^2)
  #e9 = Sys.time()
  #print(e9-s9)
  #SSE <- SSE*100
  #e = Sys.time()
  #print(e-s)
  return(SSE)
}

############################################################
SSE2 <- function(theta,
                indRainDate,
                paramGR4J,
                virObsFlow){
  #s = Sys.time()
  #Passing element in theta to WGEN parameter
  #Occurence model parameters
  #s1 = Sys.time()
  paramMC <- data.frame(matrix(NA,12,2))
  paramMC[,1] <- theta[1:12]; paramMC[,2] <- theta[13:24]
  #e1 = Sys.time()
  #print(e1-s1)
  #Amount model parameters
  #s2 = Sys.time()
  paramAmount <- data.frame(matrix(NA,12,2))
  paramAmount[,1] <- theta[25:36]; paramAmount[,2] <- theta[37:48]
  #e2 = Sys.time()
  #print(e2-s2)
  #Generate sim rain with given parameters above (a vector)
  #s3 = Sys.time()
  simRainRep <- amountModel_V3.0(occurParam = paramMC,amountParam = paramAmount, indRainDate = indRainDate, rep = 1)
  
  #e3 = Sys.time()
  #print(e3-s3)
  #Get Rainfall replicates
  #Generate sim flow with sim rain
  #add simRain to paramGR4J options
  #s4 = Sys.time()
  paramGR4J[[3]]$Precip <- simRainRep[,1]
  #e4 = Sys.time()
  #print(e4 - s4)
  #RunGR4J model with updated sim rain
  #s5 = Sys.time()
  outputGR4J <- runGR4J(paramGR4J)
  #e5 = Sys.time()
  #print(e5-s5)
  #Get sim flow from output GR4J
  #s6 = Sys.time()
  simFlowRep <- outputGR4J$Qsim
  #e6 = Sys.time()
  #print(e6-s6)
  #Calculate Exceedance Probability for sim flow and virobs flow
  #s7 = Sys.time()
  simFDC <- getExceedProb(simFlowRep)
  #e7 = Sys.time()
  #print(e7-s7)
  
  #s8 = Sys.time()
  virObsFDC <- getExceedProb(virObsFlow)
  #e8 = Sys.time()
  #print(e8-s8)
  #Calculate the Sum of square Error
  #s9 = Sys.time()
  err <- simFDC$Flow - virObsFDC$Flow
  SSE <- sum(err^2)
  #e9 = Sys.time()
  #print(e9-s9)
  #SSE <- SSE*100
  #e = Sys.time()
  #print(e-s)
  return(SSE)
}

#------------------#
OF_1 <- quote(
  SSE_FlowDurationCurve_V1.0(
    theta = iniTheta,
    obs.data = RainDatFormat,
    paramGR4J = paramGR4J,
    virObsFlow = virObsFlow
  )
)
#------------------#
OF_2 <- quote(
  SSE_FlowDurationCurve_V2.0(
    theta = iniTheta,
    indRainDate = indRainDate,
    paramGR4J = paramGR4J,
    virObsFlow = virObsFlow
  )
)
#------------------#
virObsFDC <- getExceedProb(virObsFlow)
OF_3 <- quote(
  SSE_FlowDurationCurve_V3.0(
    theta = iniTheta,
    indRainDate = indRainDate,
    paramGR4J = paramGR4J[[1]],
    inputGR4J = paramGR4J[[3]],
    runOptionGR4J = paramGR4J[[4]],
    virObsFDC = virObsFDC$Flow
  )
)
#------------------#
OF_4 <- quote(
SSE_FlowDurationCurve_V4.0(
  theta = iniTheta,
  indRainDate = indRainDate,
  paramGR4J = paramGR4J[[1]],
  inputGR4J = paramGR4J[[3]],
  runOptionGR4J = paramGR4J[[4]],
  virObsFDC = virObsFDC$Flow
)
)
#------------------#
OF_5 <- quote(
  SSE_FlowDurationCurve_V5.0(
    theta = iniTheta,
    indRainDate = indRainDate,
    paramGR4J = paramGR4J[[1]],
    inputGR4J = paramGR4J[[3]],
    runOptionGR4J = paramGR4J[[4]],
    virObsFDC = virObsFDC$Flow
  )
)
#------------------#
OFList <- list(OF_1,OF_5)
#------------------#
res <- microbenchmark::microbenchmark(list = OFList)
#------------------#
nameList <- c("OF_1","OF_5")
#------------------#
boxplot(res, names = nameList, log = FALSE)
#------------------#
res
ggplot2::autoplot(res, log = FALSE)
#------------------#


########################################################
paramMC <- data.frame(matrix(NA,12,2))
paramMC[,1] <- iniTheta[1:12]; paramMC[,2] <- iniTheta[13:24]
#Amount model parameters
paramAmount <- data.frame(matrix(NA,12,2))
paramAmount[,1] <- iniTheta[25:36]; paramAmount[,2] <- iniTheta[37:48]
#------------------#
model_V1.0 <-
  quote(
    Amount_model(
      occur.param = paramMC,
      amount.param = paramAmount,
      obs.data = RainDatFormat,
      rep = 1
    )
  )
#------------------#
model_V2.0 <-
  quote(
    amountModel_V2.0(
      occurParam = paramMC,
      amountParam = paramAmount,
      indRainDate = indRainDate,
      rep = 1
    )
  )
#------------------#
model_V3.0 <-
  quote(
    amountModel_V3.0(
      occurParam = paramMC,
      amountParam = paramAmount,
      indRainDate = indRainDate,
      rep = 1
    )
  )
#------------------#
model_V4.0 <-
  quote(
    amountModel_V4.0(
      occurParam = paramMC,
      amountParam = paramAmount,
      indRainDate = indRainDate$i.mm,
      rep = 1
    )
  )
#------------------#
model_V1.1 <-
  quote(manualWGEN(
    paramMC = paramMC,
    paramAmount = paramAmount,
    obs.data = RainDatFormat,
    rep = 1
  ))
#
OFList <- list(model_V1.0, model_V1.1, model_V2.0, model_V3.0, model_V4.0)
#------------------#
res <- microbenchmark::microbenchmark(list = OFList)
#------------------#
nameList <- c("V1.0", "V1.1", "V2.0", "V3.0", "V4.0")
#------------------#
boxplot(res, names = nameList)
#------------------#
res
#------------------#
#
exceed_V1.0 <- quote(getExceedProb(virObsFlow))
exceed_V2.0 <- quote(getExceedProb_V2.0(virObsFlow))
exceed_V3.0 <- quote(getExceedProb_V3.0(virObsFlow))
#------------------#
OFList <- list(exceed_V1.0, exceed_V2.0, exceed_V3.0)
#------------------#
res <- microbenchmark::microbenchmark(list = OFList)
#------------------#
nameList <- c("V1.0", "V2.0", "V3.0")
#------------------#
boxplot(res, names = nameList)
#------------------#
res
#------------------#
#
head(getExceedProb_V3.0(virObsFlow))



MCmodel <- quote(MCmodel(length(RainDatFormat[RainDatFormat$month==1,2]), paramMC[1,1], paramAmount[1,2]))

runT <- microbenchmark::microbenchmark(MCmodel(length(RainDatFormat[RainDatFormat$month==1,2]), paramMC[1,1], paramAmount[1,2]))
runT



getSSE <- function(simFDC, virObsFDC) {
  err <- simFDC$Flow - virObsFDC$Flow
  SSE <- sum(err^2)
  return(SSE)
}

getSSE2 <- function(simFDC, virObsFDC){
  SSE <- sum((simFDC$Flow - virObsFDC$Flow)^2)
  return(SSE)
}
simFDC <- getExceedProb(simFlowRep[,1])

res <- microbenchmark::microbenchmark(getSSE(simFDC, virObsFDC), getSSE2(simFDC, virObsFDC))
boxplot(res)
res




simrep<-amountModel_V4.0(
  occurParam = paramMC,
  amountParam = paramAmount,
  indRainDate = indRainDate,
  rep = 100
)
simRep2 <- manualWGEN(paramMC = paramMC, paramAmount = paramAmount, obs.data = RainDatFormat, rep = 1)
simRep1 <- simrep[,1]



