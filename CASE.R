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
rep <- simRainRep[[1]] #100 column Original
repoptim <- simFlowRep_Opt[[3]] #Optim
indFlowDate <- indFlowDate <- makeObsDates(RainDat$Date_Time[paramGR4J[[2]]])

#Get CASE for original rain
CASEMonthlyTotoal<-getCASEMonthlyTotal(obs = obs, sim = rep, indRainDate = indRainDate)
CASEMonthlyWetDay <- getCASEWetDayMonthly(obs = obs, sim = rep, indDate = indRainDate)
CASEMeanDayTotal <- getCASEMeanDayTotal(obs = obs, sim = rep, indRainDate = indRainDate)
CASERainMatrixORG <- rbind(CASEMonthlyWetDay,CASEMonthlyTotoal,CASEMeanDayTotal)
write.csv(CASERainMatrixORG,file = "resultCASE/rainCASEORG_Site23.csv")

#Get CASE for optim rain
CASEMonthlyTotoal<-getCASEMonthlyTotal(obs = obs, sim =  simFlowRep_Opt[[3]], indRainDate = indRainDate)
CASEMonthlyWetDay <- getCASEWetDayMonthly(obs = obs, sim =  simFlowRep_Opt[[3]], indDate = indRainDate)
CASEMeanDayTotal <- getCASEMeanDayTotal(obs = obs, sim =  simFlowRep_Opt[[3]], indRainDate = indRainDate)
CASERainMatrixOPT <- rbind(CASEMonthlyWetDay,CASEMonthlyTotoal,CASEMeanDayTotal)
write.csv(CASERainMatrixOPT,file = "resultCASE/rainCASEOPT_Site08.csv")


#Get CASE for original flow
CASEPercFlow <- getCASEPercFlow(obsFlow = virObsFlow, simFlow = simFlowRep, indFlowDate = indFlowDate)
CASEMonthlyTotoal<-getCASEMonthlyTotalFlow(obs = virObsFlow, sim = simFlowRep, indRainDate = indFlowDate)
CASEFlowMatrixORG <- rbind(CASEPercFlow, CASEMonthlyTotoal)
write.csv(CASEFlowMatrixORG,file = "resultCASE/flowCASEORG_Site08.csv")

#Get CASE for optim flow
CASEPercFlow2 <- getCASEPercFlow(obsFlow = virObsFlow, simFlow = simFlowRep_Opt[[1]], indFlowDate = indFlowDate)
CASEMonthlyTotoal2<-getCASEMonthlyTotalFlow(obs = virObsFlow, sim = simFlowRep_Opt[[1]], indRainDate = indFlowDate)
CASEFlowMatrixOPT <- rbind(CASEPercFlow2, CASEMonthlyTotoal)
write.csv(CASEFlowMatrixOPT,file = "resultCASE/flowCASEOPT_Site08.csv")


