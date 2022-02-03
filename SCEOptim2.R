#----Using SCE Optimiser----
#REQUIRE:
###RainDat, 
###paramGR4J, 
###VirObsFlow, 
###indRainDate

#Getting original parameters of the WGEN model----
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

virObsFDC <- getExceedProb_V2.0(virObsFlow)
for (i in 1:20){
  
  #Run SCE optim
  optResult <- hydromad::SCEoptim(FUN = SSE_FlowDurationCurve_V4.0,
                                  par = iniTheta,
                                  indRainDate = indRainDate,
                                  paramGR4J = paramGR4J[[1]],
                                  inputGR4J = paramGR4J[[3]],
                                  runOptionGR4J = paramGR4J[[4]],
                                  virObsFDC = virObsFDC$flow,
                                  lower = lowerTheta,
                                  upper = upperTheta)
  
  #Collect opt result and update initial parameter for the next run
  iniTheta <- assign(paste("thetaTrial_", i, sep = ""), optResult$par)
  assign(paste("optResultTrial_", i, sep = ""), optResult)
  
}
save.image(file = "SCEOptimC.RData")

#----plotting------------------------------------------------
simFlowRep_Opt <- getSimFlowRep_Opt(theta = thetaTrialIM,
                                    paramGR4J = paramGR4J,
                                    rep = 100,
                                    obsRain = RainDat)


simRainRep_Opt <- getSimRainRep(simFlowRep_Opt[[2]])

indRainDate <- makeObsDates(RainDat[,1])

pdf(
  file = paste(
    WD,
    "/Results&Plots/SCEoptim/allParam.pdf", sep=""),
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



write.csv(thetaTrial_20, file = "theta.csv")

SSE_FlowDurationCurve_V3.0(thetaTrial_20, indRainDate, paramGR4J[[1]], paramGR4J[[3]], paramGR4J[[4]], virObsFDC$flow)


#Getting influencing month parameters of the WGEN model----
iniThetaDataframe <-
  read.csv(file = paste(WD, "/WGENmodel_parameter_eachSite/Site1_410730.csv", sep = ""))
iniTheta <- rep(0,12)
iniTheta[1:3]<-iniThetaDataframe[6:8,1] 
iniTheta[4:6]<-iniThetaDataframe[6:8,2] 
iniTheta[7:9]<-iniThetaDataframe[6:8,3] 
iniTheta[10:12]<-iniThetaDataframe[6:8,4] 

#Setting boundaries
lowerTheta <- rep(0,12)
lowerTheta[1:6] = 0; lowerTheta[7:12] = 1e-10
upperTheta <- rep(0,12)
upperTheta[1:6] = 1; upperTheta[7:12] = Inf

virObsFDC <- getExceedProb_V2.0(virObsFlow)
for (i in 21:25){
  
  #Run SCE optim
  optResult <- hydromad::SCEoptim(FUN = SSE_FlowDurationCurve_InfMonth,
                                  par = thetaTrialIM_20,
                                  indRainDate = indRainDate,
                                  paramGR4J = paramGR4J[[1]],
                                  inputGR4J = paramGR4J[[3]],
                                  runOptionGR4J = paramGR4J[[4]],
                                  virObsFDC = virObsFDC$flow,
                                  lower = lowerTheta,
                                  upper = upperTheta,
                                  paramWGEN = paramWGEN)
  
  #Collect opt result and update initial parameter for the next run
  iniTheta <- assign(paste("thetaTrialIM", i, sep = ""), optResult$par)
  assign(paste("optResultTrialIM_", i, sep = ""), optResult)
  
}



paramWGEN_IM <- paramWGEN
paramWGEN_IM[6:8,1] <- thetaTrialIM_10[1:3]
paramWGEN_IM[6:8,2] <- thetaTrialIM_10[4:6]
paramWGEN_IM[6:8,3] <- thetaTrialIM_10[7:9]
paramWGEN_IM[6:8,4] <- thetaTrialIM_10[10:12]
thetaTrialIM <- rep(0,48)
thetaTrialIM[1:12]<-paramWGEN_IM[1:12,1]; thetaTrialIM[13:24]<-paramWGEN_IM[1:12,2]
thetaTrialIM[25:36]<-paramWGEN_IM[1:12,3]; thetaTrialIM[37:48]<-paramWGEN_IM[1:12,4]
