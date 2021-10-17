##Get work directory
WD<-getwd()
##Get flow, rain and evap data
RainDat <- read.csv(paste(WD,"/Data/Catchment02/70317_SILO_Rain.csv",sep = ""))
EvapDat <- read.csv(paste(WD,"/Data/Catchment02/70317_SILO_Evap.csv",sep = ""))
FlowDat <- read.csv(paste(WD,"/Data/Catchment02/410730_HRS_FLow.csv",sep = ""))



##Fit the WGEN model and Generate some rainfall replicates
RainDatFormat <- format_TimeSeries(RainDat)
SimRainList <- getSimRain(RainDatFormat, rep = 10, mod = "gama")
Simrain <- rbind(SimRainList[[1]][2],SimRainList[[2]][2],SimRainList[[3]][2],SimRainList[[4]][2],SimRainList[[5]][2],SimRainList[[6]][2],SimRainList[[7]][2],SimRainList[[8]][2],SimRainList[[9]][2],SimRainList[[10]][2],SimRainList[[11]][2],SimRainList[[12]][2])
##---------------------------------------------------------------------#
##Fit GR4J model
inputGR4J <- makeInputGR4J(P = RainDat$Value, Q = FlowDat$Value, E = EvapDat$Value, Date = FlowDat$Date_Time, A = 130)
paramGR4J <- getParamGR4J(inputGR4J = inputGR4J)

##Simulate flow time series
outputGR4J <- runGR4J(paramGR4J)
plot(outputGR4J, Qobs = inputGR4J$Q[paramGR4J[[2]]])

