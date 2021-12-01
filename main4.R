#shuffle study

obsRainShuffle <- RainDat$Value
for (m in 1:12){
  obsRainShuffle[indRainDate$i.mm[[m]]] <- randomShuffle(obsRainShuffle[indRainDate$i.mm[[m]]],10)
}

##Set GR4J model parameter
inputGR4JShuffle <-
  makeInputGR4J(
    P = obsRainShuffle,
    Q = FlowDat$Value,
    E = EvapDat$Value,
    Date = FlowDat$Date_Time,
    A = whichSite$area
  )

#Set start and end of period
start <- FlowDat[1,1] ; end <- FlowDat[length(FlowDat[,1]),1]

#months, set warm-up period, airGR prefer more than 12 months
warmup <- 24

#make parameter list
paramGR4JShuffle <-
  getParamGR4J(
    inputGR4J = inputGR4JShuffle,
    start = start,
    end = end,
    warmup = warmup,
    parameter = "known",
    knownParam = catchmentInfo[i, c(9, 10, 11, 12)]
  )

##Run GR4J model
outputGR4JShuffle <- runGR4J(paramGR4JShuffle)

virObsFlowShuffle <- outputGR4JShuffle$Qsim

  obsFDCShuffle <- getExceedProb(virObsFlowShuffle)
  
  plotFlowDurationCurve(simFlowRep = simFlowRep, virObsFlow = virObsFlow, option = "withCILimit")
  
lines(obsFDCShuffle, col="chartreuse", lwd = 2.5)
abline(v = seq(0, 1, 0.2), h = c(0.1,0.5,2,5,20), col = "lightgray", lty = 3)
  compareAnnualMaxima(indObsDate = indFlowDate, obs = virObsFlow, simRep = simFlowRep)
obsAMShuffle <- getAnnualMaxima(virObsFlowShuffle, indObsDate = indFlowDate)

obsAMRIShuffle <- getAnnualRetInt(obsAMShuffle)
points(obsAMRIShuffle,col="chartreuse")




#Check MC model----

occParamZero <- fitMCModel(RainDat)
occParamNode1 <- fitMCModel(RainDat)
occParamZero2 <- fitMCModel(RainDat)
write.csv(RainDat,file = "RainDat.csv")

RainDatMod <- RainDat
RainDatMod$Value[RainDatMod$Value <= 0.1] <- 0
occurParamNode12 <- fitMCModel(RainDatMod)

amountparam2 <- fitAmountModel(RainDat,mod = "gama")
