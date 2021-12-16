#----Using SCE Optimiser----

iniThetaDataframe <-
  read.csv(file = paste(WD, "/WGENmodel_parameter_eachSite/Site1_410730.csv", sep = ""))
iniTheta <- rep(0,48)
iniTheta[1:12]<-iniThetaDataframe[1:12,1]; iniTheta[13:24]<-iniThetaDataframe[1:12,2]
iniTheta[25:36]<-iniThetaDataframe[1:12,3]; iniTheta[37:48]<-iniThetaDataframe[1:12,4]
#iniTheta[1:24] <- theta_TrialB01
lowerTheta <- rep(0,48)
lowerTheta[1:24] = 0; lowerTheta[25:48] = 1e-10
upperTheta <- rep(0,48)
upperTheta[1:24] = 1; upperTheta[25:48] = Inf

#iniTheta <- theta_TrialA02

opt <- hydromad::SCEoptim(FUN = SSE_FlowDurationCurve,
                          par = iniTheta,
                          obsRain = RainDat,
                          paramGR4J = paramGR4J,
                          virObsFlow = virObsFlow,
                          lower = lowerTheta,
                          upper = upperTheta)
Sys.time()



theta_TrialA03 <- opt$par
theta <- theta_TrialA03
optTrialA03 <- opt
##Collecting results from trial
trial <- list(theta_Trial01, theta_Trial02, theta_Trial03,theta_Trial04,theta_Trial05,theta_Trial06,theta_Trial07,theta_Trial08)

for (t in 6:8){
  theta <- trial[[t]]
  
  simFlowRep_Opt <- getSimFlowRep_Opt(theta = thetaMC,
                                      paramGR4J = paramGR4J,
                                      rep = 100,
                                      obsRain = RainDat)
  
  
  simRainRep_Opt <- getSimRainRep(simFlowRep_Opt[[2]])
  indRainDate <- makeObsDates(RainDat$Date_Time[paramGR4J[[2]]]) ##Get date index
  
  pdf(
    file = paste(
      WD,
      "/Results&Plots/SCEoptim/allParam_",
      "Site",
      s,
      "_",
      whichSite$station,
      "_",
      "Trial",
      t,
      ".pdf",
      sep = ""
    ),
    width = 8.25,
    height = 11.75
  )
  
  par(mfrow=c(2,2))
  
  compareAnnualMaxima(indObsDate = indRainDate, obs = RainDat$Value, simRep = simRainRep_Opt)
  
  siloStation <-
    paste(
      " Daily Annual Maxima (Rainfall)",
      "\n",
      " ID:",
      siloInfo$station,
      "Name:",
      as.character(siloInfo$name),
      "\n",
      "Juradiction:",
      siloInfo$state,
      "Lat:",
      siloInfo$latitude,
      "Long:",
      siloInfo$longitude,
      "Elevation:",
      siloInfo$elevation,
      sep = " "
    )
  
  mtext(siloStation,
        3,
        0,
        cex = 0.5,
        adj = 0,
        padj = -0.3)
  
  plotFlowDurationCurve(simFlowRep = simFlowRep_Opt[[1]],
                        virObsFlow = virObsFlow,
                        option = "withCILimit")
  
  outlet <-
    paste(
      " Flow Duration Curve",
      "\n",
      " ID:",
      whichSite[1],
      "Juradiction:",
      whichSite[4],
      "Area:",
      whichSite[5],
      "km^2",
      "\n",
      "Lat:",
      whichSite[2],
      "Long:",
      whichSite[3],
      "From:",
      start,
      "To:",
      end,
      sep = " "
    )
  
  mtext(outlet,
        3,
        0,
        cex = .5,
        adj = 0,
        padj = -.3)
  indFlowDate <- makeObsDates(RainDat$Date_Time[paramGR4J[[2]]])
  compareAnnualMaxima(indObsDate = indFlowDate,
                      obs = virObsFlow,
                      simRep = simFlowRep_Opt[[1]])
  
  mtext("Daily annual maxima (Flow)")
  RainDatFormat <- format_TimeSeries(RainDat)
  wetday_monthlystats_plot(RainDatFormat, simFlowRep_Opt[[2]], type = "boxplot")
  monthlytotal_stats_plot(RainDatFormat, simFlowRep_Opt[[2]], type = "boxplot")
  compareMean3dayTotal(obs = RainDat$Value, sim = simFlowRep_Opt[[3]], indObsDate = indRainDate)
  compareMean5dayTotal(obs = RainDat$Value, sim = simFlowRep_Opt[[3]], indObsDate = indRainDate)
  #Dry Spell Plot
  month <-
    c("Jan",
      "Feb",
      "Mar",
      "Apr",
      "May",
      "Jun",
      "Jul",
      "Aug",
      "Sep",
      "Oct",
      "Nov",
      "Dec")
  
  maxSpell <- 10
  
  for (m in 1:12){
    compareDrySpell(
      monthlyObsRain = RainDat$Value[indRainDate$i.mm[[m]]],
      monthlySimRain = simRainRep_Opt[indRainDate$i.mm[[m]], ],
      maxSpell = maxSpell
    )
    mtext(paste("Dry Spell: ", paste(month[m]), sep=""))
    
    #Wet Spell Plot
    compareWetSpell(
      monthlyObsRain = RainDat$Value[indRainDate$i.mm[[m]]],
      monthlySimRain = simRainRep_Opt[indRainDate$i.mm[[m]], ],
      maxSpell = maxSpell
    )
    mtext(paste("Wet Spell: ", paste(month[m]), sep=""))
  }
  
  dev.off()
}



paramocc <- data.frame(matrix(NA,12,2))
paramocc[,1] <- theta_TrialA09[1:12]; paramocc[,2] <- theta_TrialA09[13:24]
paramamo <- data.frame(matrix(NA,12,2))
paramamo[,1] <- theta_TrialA09[25:36]; paramamo[,2] <- theta_TrialA09[37:48]

write.csv(paramocc,"occur.csv")
write.csv(paramamo,"amount.csv")


#SCEOptim occurence parameter only
paramAmount <- rep(0,24)
paramAmount[1:12] <- iniThetaDataframe[1:12,3]; paramAmount[13:24]<-iniThetaDataframe[1:12,4]

iniThetaDataframe <-
  read.csv(file = paste(WD, "/WGENmodel_parameter_eachSite/Site1_410730.csv", sep = ""))
iniTheta <- rep(0,24)
iniTheta[1:12]<-iniThetaDataframe[1:12,1]; iniTheta[13:24]<-iniThetaDataframe[1:12,2]
#iniTheta[25:36]<-iniThetaDataframe[1:12,3]; iniTheta[37:48]<-iniThetaDataframe[1:12,4]

lowerTheta <- rep(0,24)
lowerTheta[1:24] = 0;
upperTheta <- rep(0,24)
upperTheta[1:24] = 1

SSE_FlowDurationCurveOcc <- function(theta,
                                  obsRain,
                                  paramGR4J,
                                  virObsFlow,
                                  occ){
  #Passing element in theta to WGEN parameter
  #Occurence model parameters
  paramMC <- data.frame(matrix(NA,12,2))
  paramMC[,1] <- theta[1:12]; paramMC[,2] <- theta[13:24]
  #Amount model parameters
  paramAmount <- data.frame(matrix(NA,12,2))
  paramAmount[,1] <- occ[1:12]; paramAmount[,2] <- occ[13:24]
  
  #Generate sim rain with given parameters above (a vector)
  simRainRep <- manualWGEN(paramMC, paramAmount, obsRain, rep = 1) #Get Rainfall replicates
  #Generate sim flow with sim rain
  #add simRain to paramGR4J options
  paramGR4J[[3]]$Precip <- simRainRep
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

iniTheta <- thetaMC_Trial14
opt <- hydromad::SCEoptim(FUN = SSE_FlowDurationCurveOcc,
                          par = iniTheta,
                          obsRain = RainDat,
                          paramGR4J = paramGR4J,
                          virObsFlow = virObsFlow,
                          occ = paramAmount,
                          lower = lowerTheta,
                          upper = upperTheta)
Sys.time()
#------------------------------------------------------#
save.image(file = "SCEOptimMC.RData")

#--------------------------------------------------

sse <- rep (0,100)
for (i in 1:100){
  sse[i] <- SSE_FlowDurationCurveOcc(theta = thetaMC_Trial19,
                                  obsRain = RainDat,
                                  paramGR4J = paramGR4J,
                                  virObsFlow = virObsFlow,
                                  occ = paramAmount)
}
sse <- mean(sse)

sse <- rep (0,100)
for (i in 1:100){
  sse[i] <- SSE_FlowDurationCurveAmo(theta = thetaAmo_Trial10,
                                     obsRain = RainDat,
                                     paramGR4J = paramGR4J,
                                     virObsFlow = virObsFlow,
                                     amo = paramMC)
}
sse <- mean(sse)

sse <- rep (0,100)
for (i in 1:100){
  sse[i] <- SSE_FlowDurationCurve(theta = iniTheta,
                                     obsRain = RainDat,
                                     paramGR4J = paramGR4J,
                                     virObsFlow = virObsFlow
                                     )
}
sse <- mean(sse)


thetaMC <- rep(0,48)
thetaMC[1:24] <- thetaMC_Trial19
thetaMC[25:48] <- paramAmount

#--------------------------------------------------
#SCEOptim Amount parameter only with optimized MC param
paramMC <- thetaMC_Trial16

iniThetaDataframe <-
  read.csv(file = paste(WD, "/WGENmodel_parameter_eachSite/Site1_410730.csv", sep = ""))
iniTheta <- rep(0,24)
iniTheta[1:12]<-iniThetaDataframe[1:12,3]; iniTheta[13:24]<-iniThetaDataframe[1:12,4]

lowerTheta <- rep(0,24)
lowerTheta[1:24] = 1e-10;
upperTheta <- rep(0,24)
upperTheta[1:24] = Inf

SSE_FlowDurationCurveAmo <- function(theta,
                                     obsRain,
                                     paramGR4J,
                                     virObsFlow,
                                     amo){
  #Passing element in theta to WGEN parameter
  #Occurence model parameters
  paramMC <- data.frame(matrix(NA,12,2))
  paramMC[,1] <- amo[1:12]; paramMC[,2] <- amo[13:24]
  #Amount model parameters
  paramAmount <- data.frame(matrix(NA,12,2))
  paramAmount[,1] <- theta[1:12]; paramAmount[,2] <- theta[13:24]
  
  #Generate sim rain with given parameters above (a vector)
  simRainRep <- manualWGEN(paramMC, paramAmount, obsRain, rep = 1) #Get Rainfall replicates
  #Generate sim flow with sim rain
  #add simRain to paramGR4J options
  paramGR4J[[3]]$Precip <- simRainRep
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

opt <- hydromad::SCEoptim(FUN = SSE_FlowDurationCurveAmo,
                          par = iniTheta,
                          obsRain = RainDat,
                          paramGR4J = paramGR4J,
                          virObsFlow = virObsFlow,
                          amo = paramMC,
                          lower = lowerTheta,
                          upper = upperTheta)
Sys.time()
#------------------------------------------------------#
thetaAmo_Trial01 <- opt$par
optAmo_Trial01 <- opt
iniTheta <- thetaAmo_Trial01
opt <- hydromad::SCEoptim(FUN = SSE_FlowDurationCurveAmo,
                          par = iniTheta,
                          obsRain = RainDat,
                          paramGR4J = paramGR4J,
                          virObsFlow = virObsFlow,
                          amo = paramMC,
                          lower = lowerTheta,
                          upper = upperTheta)
Sys.time()
#------------------------------------------------------#
thetaAmo_Trial02 <- opt$par
optAmo_Trial02 <- opt
iniTheta <- thetaAmo_Trial02
opt <- hydromad::SCEoptim(FUN = SSE_FlowDurationCurveAmo,
                          par = iniTheta,
                          obsRain = RainDat,
                          paramGR4J = paramGR4J,
                          virObsFlow = virObsFlow,
                          amo = paramMC,
                          lower = lowerTheta,
                          upper = upperTheta)
Sys.time()
#------------------------------------------------------#
thetaAmo_Trial03 <- opt$par
optAmo_Trial03 <- opt
iniTheta <- thetaAmo_Trial03
opt <- hydromad::SCEoptim(FUN = SSE_FlowDurationCurveAmo,
                          par = iniTheta,
                          obsRain = RainDat,
                          paramGR4J = paramGR4J,
                          virObsFlow = virObsFlow,
                          amo = paramMC,
                          lower = lowerTheta,
                          upper = upperTheta)
Sys.time()
#------------------------------------------------------#
thetaAmo_Trial04 <- opt$par
optAmo_Trial04 <- opt
iniTheta <- thetaAmo_Trial04
opt <- hydromad::SCEoptim(FUN = SSE_FlowDurationCurveAmo,
                          par = iniTheta,
                          obsRain = RainDat,
                          paramGR4J = paramGR4J,
                          virObsFlow = virObsFlow,
                          amo = paramMC,
                          lower = lowerTheta,
                          upper = upperTheta)
Sys.time()
#------------------------------------------------------#
thetaAmo_Trial05 <- opt$par
optAmo_Trial05 <- opt
iniTheta <- thetaAmo_Trial05
opt <- hydromad::SCEoptim(FUN = SSE_FlowDurationCurveAmo,
                          par = iniTheta,
                          obsRain = RainDat,
                          paramGR4J = paramGR4J,
                          virObsFlow = virObsFlow,
                          amo = paramMC,
                          lower = lowerTheta,
                          upper = upperTheta)
Sys.time()
#------------------------------------------------------#
thetaAmo_Trial06 <- opt$par
optAmo_Trial06 <- opt
iniTheta <- thetaAmo_Trial06
opt <- hydromad::SCEoptim(FUN = SSE_FlowDurationCurveAmo,
                          par = iniTheta,
                          obsRain = RainDat,
                          paramGR4J = paramGR4J,
                          virObsFlow = virObsFlow,
                          amo = paramMC,
                          lower = lowerTheta,
                          upper = upperTheta)
Sys.time()
#------------------------------------------------------#
thetaAmo_Trial07 <- opt$par
optAmo_Trial07 <- opt
iniTheta <- thetaAmo_Trial07
opt <- hydromad::SCEoptim(FUN = SSE_FlowDurationCurveAmo,
                          par = iniTheta,
                          obsRain = RainDat,
                          paramGR4J = paramGR4J,
                          virObsFlow = virObsFlow,
                          amo = paramMC,
                          lower = lowerTheta,
                          upper = upperTheta)
Sys.time()
#------------------------------------------------------#
thetaAmo_Trial08 <- opt$par
optAmo_Trial08 <- opt
iniTheta <- thetaAmo_Trial08
opt <- hydromad::SCEoptim(FUN = SSE_FlowDurationCurveAmo,
                          par = iniTheta,
                          obsRain = RainDat,
                          paramGR4J = paramGR4J,
                          virObsFlow = virObsFlow,
                          amo = paramMC,
                          lower = lowerTheta,
                          upper = upperTheta)
Sys.time()
#------------------------------------------------------#
thetaAmo_Trial09 <- opt$par
optAmo_Trial09 <- opt
iniTheta <- thetaAmo_Trial09
opt <- hydromad::SCEoptim(FUN = SSE_FlowDurationCurveAmo,
                          par = iniTheta,
                          obsRain = RainDat,
                          paramGR4J = paramGR4J,
                          virObsFlow = virObsFlow,
                          amo = paramMC,
                          lower = lowerTheta,
                          upper = upperTheta)
Sys.time()
#------------------------------------------------------#
thetaAmo_Trial10 <- opt$par
optAmo_Trial10 <- opt
save.image(file = "SCEOptimAmo.RData")
thetaNew <- rep(0,48)
thetaNew[1:24] <- paramMC
thetaNew[25:48] <- thetaAmo_Trial10
