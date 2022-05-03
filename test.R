i=11
simFDC <- getExceedProb_V2.0(simFlowRep[indFlowDate$i.mm[[i]],1])
virObsFDC <- getExceedProb_V2.0(virObsFlow[indFlowDate$i.mm[[i]]])



err1 <- abs(((simFDC$flow[which(simFDC$Prob<=0.5)] - virObsFDC$flow[which(virObsFDC$Prob<=0.5)]))/virObsFDC$flow[which(virObsFDC$Prob<=0.5)])

err2 <- abs(((simFDC$flow[which(simFDC$Prob>=0.5)] - virObsFDC$flow[which(virObsFDC$Prob>=0.5)]))/virObsFDC$flow[which(virObsFDC$Prob>=0.5)])



SRE1 <- sum(err1)
SRE2 <- sum(err2)
par(mfrow = c(3,4))
for (i in (1:12)){
  plotFlowDurationCurve(simFlowRep[indFlowDate$i.mm[[i]],],virObsFlow[indFlowDate$i.mm[[i]]],"withCILimit")
}



layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plotFlowDurationCurve(simFlowRep = simFlowRep_Opt[[1]],
                      virObsFlow = virObsFlow,
                      option = "withCILimit")
plotFlowPercentilesV2.0(obs = virObsFlow, sim = simFlowRep, indFlowDate = indFlowDate, optimSim = simFlowRep_Opt[[1]], mod="2", percentile = c(0.05, 0.5, 0.95))
plotFlowPercentiles(obs = virObsFlow, sim = simFlowRep, indFlowDate = indFlowDate, percentile = c(0.05, 0.95))

par(mfrow=c(1,2))
plotFlowDurationCurve(simFlowRep = simFlowRep_Opt[[1]],
                      virObsFlow = virObsFlow,
                      option = "withCILimit")
plotFlowPercentilesV2.0(obs = virObsFlow, sim = simFlowRep, indFlowDate = indFlowDate, optimSim = simFlowRep_Opt[[1]], mod="2", percentile = c(0.05, 0.5, 0.95))


write.csv(simFlowRep, file = "fixedSimFlowCatch08.csv")
write.csv(simRainRep[[1]],file = "fixedSimRainCatch08.csv")
fixedSimFlowCatch08 <- read.csv(file = "fixedSimFlowCatch08.csv")
fixedSimRainCatch08 <- read.csv(file = "fixedSimRainCatch08.csv")


write.csv(simFlowRep, file = "fixedSimFlowCatch23.csv",row.names = F)
write.csv(simRainRep[[1]],file = "fixedSimRainCatch23.csv",row.names = F)
fixedSimFlowCatch23 <- read.csv(file = "fixedSimFlowCatch23.csv")
fixedSimRainCatch23 <- read.csv(file = "fixedSimRainCatch23.csv")
plotFlowPercentiles(obs = virObsFlow, sim = fixedSimFlowCatch23, indFlowDate = indFlowDate, percentile = c(0.05, 0.95))
CASEPercFlow <- getCASEPercFlow(obsFlow = virObsFlow, simFlow = fixedSimFlowCatch23, indFlowDate = indFlowDate)

CASEMonthlyTotal<-getCASEMonthlyTotal(obs = virObsFlow, sim = fixedSimFlowCatch08, indRainDate = indFlowDate)

simFlowRep_Opt <- getSimFlowRep_Opt(theta = thetaTrialSepmonth,
                                    paramGR4J = paramGR4J,
                                    rep = 100,
                                    obsRain = RainDat,seed)
par(mfrow=c(1,3))
plotFlowPercentilesV2.0(obs = virObsFlow, sim = simFlowRep, indFlowDate = indFlowDate, optimSim = simFlowRep_Opt[[1]], mod="2", percentile = c(0.05, 0.5, 0.95))

RainCASETotal <- getCASEMonthlyTotal(obs = RainDat[,2], sim = simFlowRep_Opt[[3]], indRainDate = indRainDate)



layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))

plotFlowDurationCurveV2.0(simFlowRep = fixedSimFlowCatch08,simFlowRepOpt = simFlowRep_Opt[[1]] ,
                          virObsFlow = virObsFlow)
plotFlowPercentilesV2.0(obs = virObsFlow, sim = fixedSimFlowCatch08, indFlowDate = indFlowDate, optimSim = simFlowRep_Opt[[1]], mod="2", percentile = c(0.05, 0.5, 0.95))
