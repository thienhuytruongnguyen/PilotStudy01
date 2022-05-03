for (i in c(16,23)){
  
  s = i
  
  outSim <- mainSimulation(i,s)
  
  #Get CASE for original rain
  #CASEMonthlyTotoal<-getCASEMonthlyTotal(obs = outSim[[1]], sim = outSim[[2]], indRainDate = outSim[[5]])
  #write.csv(CASEMonthlyTotoal,file = paste("resultCASE/rainCASEORG_Site2", i, ".csv")) 
  
  #Get CASE for original flow
  #CASEPercFlow <- getCASEPercFlow(obsFlow = outSim[[3]], simFlow = outSim[[4]], indFlowDate = outSim[[6]])
  #write.csv(CASEPercFlow,file = paste("resultCASE/flowCASEORG_Site", i, ".csv"))
  CASETotalFlow <- getCASEMonthlyTotal(obs = outSim[[3]], sim = outSim[[4]], indRainDate = outSim[[6]])
  write.csv(CASETotalFlow,file = paste("resultCASE/flowCASEORGTotal_Site", i, ".csv"))
}
