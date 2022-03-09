#------------------Runoff Model----------------------------------
optimParamWGEN <- paramWGEN
  
for (k in 1:5){
  for (i in c(2,seq(4,10))){
    ##Set GR4J model parameter
    inputGR4J_SepMonth <-
      makeInputGR4J(
        P = RainDat$Value[indRainDate$i.mm[[i]]],
        Q = FlowDat$Value[indRainDate$i.mm[[i]]],
        E = EvapDat$Value[indRainDate$i.mm[[i]]],
        Date = FlowDat$Date_Time[indRainDate$i.mm[[i]]],
        A = whichSite$area
      )
    
    #Set start and end of period
    start <- FlowDat[indRainDate$i.mm[[i]],1][1] ; end <- FlowDat[indRainDate$i.mm[[i]],1][length(indRainDate$i.mm[[i]])]
    
    
    #make parameter list
    paramGR4J_SepMonth <-
      getParamGR4J_SepMonth(
        inputGR4J = inputGR4J_SepMonth,
        start = start,
        end = end,
        parameter = "known",
        knownParam = catchmentInfo[1, c(9, 10, 11, 12)]
      )
    
    ##Run GR4J model
    outputGR4J_SepMonth <- runGR4J(paramGR4J_SepMonth)
    virObsFlowSM <- outputGR4J_SepMonth$Qsim
    virObsFDCSM <- getExceedProb_V2.0(virObsFlowSM)
    
    #Initiate parameter to optimise
    iniTheta <- rep(0,4) 
    iniTheta[1] <- optimParamWGEN[i,1]; iniTheta[2] <-optimParamWGEN[i,2];
    iniTheta[3] <- optimParamWGEN[i,3]; iniTheta[4] <- optimParamWGEN[i,4]
    lowerTheta <- rep(0,4)
    lowerTheta[1:2] = 0; lowerTheta[3:4] = 1e-10
    upperTheta <- rep(0,4)
    upperTheta[1:2] = 1; upperTheta[3:4] = Inf
    
    #Run SCE optim
    optResult <- hydromad::SCEoptim(FUN = NSE_WeightedFDC_SingleMonthRE,
                                    par = iniTheta,
                                    obsRain = RainDat[indRainDate$i.mm[[i]],2],
                                    paramGR4J = paramGR4J_SepMonth[[1]],
                                    inputGR4J = paramGR4J_SepMonth[[3]],
                                    runOptionGR4J = paramGR4J_SepMonth[[4]],
                                    virObsFlow = virObsFlowSM,
                                    lower = lowerTheta,
                                    upper = upperTheta)
    
    optimParamWGEN[i,] <- optResult$par
    remove(paramGR4J_SepMonth, inputGR4J_SepMonth)
  }
}


thetaTrialSepmonth <- rep(0,48)
thetaTrialSepmonth[1:12]<-optimParamWGEN[1:12,1]; thetaTrialSepmonth[13:24]<-optimParamWGEN[1:12,2]
thetaTrialSepmonth[25:36]<-optimParamWGEN[1:12,3]; thetaTrialSepmonth[37:48]<-optimParamWGEN[1:12,4]

#----plotting------------------------------------------------
simFlowRep_Opt <- getSimFlowRep_Opt(theta = thetaTrialSepmonth,
                                    paramGR4J = paramGR4J,
                                    rep = 100,
                                    obsRain = RainDat,seed)


simRainRep_Opt <- getSimRainRep(simFlowRep_Opt[[2]])

indRainDate <- makeObsDates(RainDat[,1])

pdf(
  file = paste(
    WD,
    "/Results&Plots/SCEoptim/OF2.pdf", sep=""),
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
plotFlowDurationCurve(simFlowRep = simFlowRep,
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
plotFlowPercentilesV2.0(obs = virObsFlow, sim = simFlowRep, indFlowDate = indFlowDate, optimSim = simFlowRep_Opt[[1]], mod="2", percentile = c(0.05, 0.5, 0.95))
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
    maxSpell = maxSpell)
}

dev.off()




#----Separate month----

Theta <- rep(0,4)
Theta[1] <- optimParamWGEN[i,1]; Theta[2] <-optimParamWGEN[i,2];
Theta[3] <- optimParamWGEN[i,3]; Theta[4] <- optimParamWGEN[i,4]

simRainSM <- WGEN_SingleMonth(occurParam = Theta[1:2], amountParam = Theta[3:4], obsRain = RainDat[indRainDate$i.mm[[4]],2],rep = 100)

#Get sim flow from output GR4J
simFlowRepSM <- getSimFlowRep(simRainRep = simRainSM, paramGR4J = paramGR4J_SepMonth)
plotFlowDurationCurve(simFlowRep = simFlowRepSM, virObsFlow = virObsFlowSM, option = "withCILimit")

percentile95(virObsFlowSM)
