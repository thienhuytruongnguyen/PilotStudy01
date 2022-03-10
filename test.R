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
plotFlowPercentiles(obs = virObsFlow, sim = simFlowRep, indFlowDate = indFlowDate, percentile = c(0.05, 0.5, 0.95))

par(mfrow=c(1,1))
plotFlowDurationCurve(simFlowRep = simFlowRep_Opt[[1]],
                      virObsFlow = virObsFlow,
                      option = "withCILimit")
plotFlowPercentilesV2.0(obs = virObsFlow, sim = simFlowRep, indFlowDate = indFlowDate, optimSim = simFlowRep_Opt[[1]], mod="2", percentile = c(0.05, 0.5, 0.95))
