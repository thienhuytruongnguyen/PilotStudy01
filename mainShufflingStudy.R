#Random shuffle----
simRainListShuffle <- SimRainList
for (m in 1:12){
  for (r in 1:100){
    simRainListShuffle[[m]][r,] <- randomShuffle(simRainListShuffle[[m]][r,],10)
  }
}

simRainRepShuffle <- getSimRainRep(simRainListShuffle)

simFlowRepShuffle <-
  getSimFlowRep(simRainRep = simRainRepShuffle, paramGR4J = paramGR4J)

plotFlowDurationCurve(simFlowRep = simFlowRep, virObsFlow = virObsFlow, option = "withCILimit")
plotFlowDurationCurve(simFlowRep = simFlowRepShuffle, virObsFlow = virObsFlow, option = "withCILimit")
compareAnnualMaxima(indObsDate = indFlowDate, obs = virObsFlow, simRep = simFlowRep)
compareAnnualMaxima(indObsDate = indFlowDate, obs = virObsFlow, simRep = simFlowRepShuffle)
#wet day amount stats
wetday_monthlystats_plot(RainDatFormat, simRainListShuffle, type = "boxplot")

#Monthly total stats
monthlytotal_stats_plot(RainDatFormat, SimRainList, type = "boxplot")

#mean 3-5 day total monthly
compareMean3dayTotal(obs = RainDat$Value, sim = simRainRep, indObsDate = indRainDate)
compareMean5dayTotal(obs = RainDat$Value, sim = simRainRep, indObsDate = indRainDate)

