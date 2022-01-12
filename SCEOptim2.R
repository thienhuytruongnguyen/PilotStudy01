#----Using SCE Optimiser----
#REQUIRE:
###RainDat, 
###paramGR4J, 
###VirObsFlow, 
###indRainDate

#Getting original parameters of the WGEN model
iniThetaDataframe <-
  read.csv(file = paste(WD, "/WGENmodel_parameter_eachSite/Site1_410730.csv", sep = ""))
iniTheta <- rep(0,48)
iniTheta[1:12]<-iniThetaDataframe[1:12,1]; iniTheta[13:24]<-iniThetaDataframe[1:12,2]
iniTheta[25:36]<-iniThetaDataframe[1:12,3]; iniTheta[37:48]<-iniThetaDataframe[1:12,4]

#Setting boundaries
lowerTheta <- rep(0,48)
lowerTheta[1:24] = 0; lowerTheta[25:48] = 1e-10
upperTheta <- rep(0,48)
upperTheta[1:24] = 1; upperTheta[25:48] = Inf


for (i in 6:12){
  
  #Run SCE optim
  optResult <- hydromad::SCEoptim(FUN = SSE_FlowDurationCurve,
                                  par = iniTheta,
                                  indRainDate = indRainDate,
                                  paramGR4J = paramGR4J,
                                  virObsFlow = virObsFlow,
                                  lower = lowerTheta,
                                  upper = upperTheta)
  
  #Collect opt result and update initial parameter for the next run
  iniTheta <- assign(paste("thetaTrial_", i, sep = ""), optResult$par)
  assign(paste("optResultTrial_", i, sep = ""), optResult)
  
}
save.image(file = "SCEOptimC.RData")
