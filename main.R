##Get work directory
runStart <- Sys.time()
WD<-getwd()
  
    ##Get site list
    flowSiteList <- read.table("siteList.txt",header=T,sep = "",dec = ".")
    weatherSiteList <- getWeatherSiteList()
    whichSite <- flowSiteList[11,]##1st row of flow site list
    
    ##Find Nearest Weather station to Site
    nearestStation <- getNearestWeatherStation(flowSiteList, weatherSiteList, whichSite = whichSite$ID)
    
    ##Getting rain and PET data from URL
    siloData <- getWeatherData(whichStation, start = whichSite$start, finish = whichSite$finish)
    
    ##Get flow, rain and evap data
    RainDat <- siloData[,c(2,3)]; colnames(RainDat) = c("Date_Time","Value");
    EvapDat <- siloData[,c(2,5)]; colnames(EvapDat) = c("Date_Time","Value")
    FlowDat <- read.csv(paste(WD,"/Data/Flow/",whichSite$station,".csv",sep = ""))
    
    ##Fit GR4J model
    inputGR4J <- makeInputGR4J(P = RainDat$Value, Q = FlowDat$Value, E = EvapDat$Value, Date = FlowDat$Date_Time, A = whichSite$area)
    start <- FlowDat[1,1] ; end <- FlowDat[length(FlowDat[,1]),1] #Set start and end of period
    warmup <- 24 #months, set warm-up period, airGR prefer more than 12 months
    paramGR4J <- getParamGR4J(inputGR4J = inputGR4J, start = start, end = end, warmup = warmup)
    
    ##Run GR4J model
    outputGR4J <- runGR4J(paramGR4J)
    
    #----------------------------------------------------#
    
    ##Fit the WGEN model and Generate some rainfall replicates
    RainDatFormat <- format_TimeSeries(RainDat)
    rep = 30 #set number of replicates
    SimRainList <- getSimRain(RainDatFormat, rep = rep, mod = "gama")
    
    ##Get SimRain Rep
    simRainRep <- getSimRainRep(SimRainList)
    
    ##Get SimFlow Rep
    simFlowRep <- getSimFlowRep(simRainRep = simRainRep, paramGR4J = paramGR4J)
    
    #Get virtual observed flow
    virObsFlow <- outputGR4J$Qsim
    
    #Plot FDC
    plotFlowDurationCurve(simFlowRep = simFlowRep, virObsFlow = virObsFlow, option = "withCILimit")
    
  
runEnd <- Sys.time()

runEnd - runStart

# #Get exceedance probability for sim flow reps
# simExceedProbRep <- getExceedProbRep(simFlowRep = simFlowRep)
# 
# #Get exceedance probability for virtual observed flow
# virObsExceedProb <- getExceedProb(flow = virObsFlow)
# 
# ##Plot

# #Plot simulated vs observed rainfall stats
# wetday_monthlystats_plot(RainDatFormat,SimRainList, type = "boxplot") #plot monthly wet day stats, type of plot can be "boxplot" or "errorbar"
# monthlytotal_stats_plot(RainDatFormat,SimRainList, type = "boxplot") #plot monthly total stats, type of plot can be "boxplot" or "errorbar"
# 
# #plot simulated vs observed flow stats
plot(outputGR4J, Qobs = inputGR4J$Q[paramGR4J[[2]]])+#plot comparison obs vs simflow_obsRain
mtext("Obseved Flow vs Simulated Flow with Observed Rain", side = 4, line = 0.5)

paramGR4J_simRain <- paramGR4J

paramGR4J_simRain[[3]]$Precip <- simRainRep[,1]

outputGR4J_simRain <- runGR4J(paramGR4J_simRain)

plot(outputGR4J_simRain, Qobs = inputGR4J$Q[paramGR4J_simRain[[2]]])#plot comparison obs vs simflow_simRain
mtext("Obseved Flow vs Simulated Flow with Simulated Rain", side = 4, line = 0.5)

#Plot FDC
plotFlowDurationCurve(simFlowRep = simFlowRep, virObsFlow = virObsFlow, option = "RepVSVirObs")




