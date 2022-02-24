## Calculating good, fair, poor degree of replicative data
#Good: <10% outside 90PL or <= 90% inside 90PL
#Fair: >10% outside 90PL but within 99.7PL or abs relative different between obs and sim mean is less than 5%
#Bad: otherwise

## >> Things neeed to calculate:
#1. 90 limits of the sim
#2. 99.7 limits of the sim
#3. absolute relative different between obs and sim

#Testing with mean monthly wet day amount

obs <- RainDat[,2] # 1 column
rep <- simRainRep[[1]] #100 column

# CASEMonthlyTotoal<-getCASEMonthlyTotal(obs = obs, sim = rep, indRainDate = indRainDate)
# CASEMonthlyWetDay <- getCASEWetDayMonthly(obs = obs, sim = rep, indDate = indRainDate)
CASEMeanDayTotal <- getCASEMeanDayTotal(obs = obs, sim = rep, indRainDate = indRainDate)
# write.csv(CASEMonthlyTotoal,file = "MT.csv")
# write.csv(CASEMonthlyWetDay,file = "WD.csv")
write.csv(CASEMeanDayTotal, file = "MD.csv")

indFlowDate <- indFlowDate <- makeObsDates(RainDat$Date_Time[paramGR4J[[2]]])
CASEPercFlow <- getCASEPercFlow(obsFlow = virObsFDC, simFlow = simFlowRep, indFlowDate = indFlowDate)
write.csv(CASEPercFlow,file = "PF.csv")

