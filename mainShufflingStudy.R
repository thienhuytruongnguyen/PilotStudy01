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






pdf(
  file = paste(
    WD,
    "/Results&Plots/RainfallStats_FlowDurationCurve_EachSite/",
    "EvapShuffling.pdf",
    sep = ""
  ),
  width = 8.25,
  height = 11.75
)

par(mfrow=c(2,2))

s=0

# Run for each of 25 sites----
for (i in 1:25) {
  
  #Updating site
  s=s+1
  
  
  #------------------Get Data----------------------------------
  
  ##Get site list
  whichSite <- flowSiteList[i,]##1st row of flow site list
  
  ##Get Nearest Weather station to Site
  siloInfo <-
    weatherSiteList[which(weatherSiteList$station == catchmentInfo$NearestSilo[i]), ]
  
  ##Getting rain and PET data from URL
  siloData <-
    getWeatherData(
      nearestStation = siloInfo$station,
      start = whichSite$start,
      finish = whichSite$finish
    )
  
  ##Get flow, rain and evap data
  RainDat <- siloData[,c(2,3)]; colnames(RainDat) = c("Date_Time","Value")
  EvapDat <- siloData[,c(2,5)]; colnames(EvapDat) = c("Date_Time","Value")
  FlowDat <- read.csv(paste(WD,"/Data/Flow/",whichSite$station,".csv",sep = ""))
  
  #------------------Runoff Model----------------------------------
  
  ##Set GR4J model parameter
  inputGR4J <-
    makeInputGR4J(
      P = RainDat$Value,
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
  paramGR4J <-
    getParamGR4J(
      inputGR4J = inputGR4J,
      start = start,
      end = end,
      warmup = warmup,
      parameter = "known",
      knownParam = catchmentInfo[i, c(9, 10, 11, 12)]
    )
  
  ##Run GR4J model
  outputGR4J <- runGR4J(paramGR4J)
  effiReport <-
    getEffiReport(
      outputGR4J = outputGR4J,
      inputGR4J = inputGR4J,
      start = start,
      end = end,
      warmup = warmup
    )
  
  
  ##Fit the WGEN model and Generate some rainfall replicates
  indRainDate <- makeObsDates(RainDat[,1]) #get index rain day
  
  virObsFlow <- outputGR4J$Qsim
  
#Make shuffled Evap dataframe
EvapShuffle <- data.frame(matrix(NA,nrow = nrow(EvapDat), ncol = 100))

for (m in 1:12){
  for (y in 1:length(indRainDate$i.ym)){
    for (r in 1:100){
      EvapShuffle[indRainDate$i.ym[[y]][[m]],r] <- randomShuffle(EvapDat$Value[indRainDate$i.ym[[y]][[m]]],10)
    }
   }
  }

#Make shuffled virobs flow
virObsFlowShuffle <- data.frame(matrix(NA,length(virObsFlow),100))
#------------------Runoff Model----------------------------------
for (r in 1:100){
  ##Set GR4J model parameter
  inputGR4JShuffle <-
    makeInputGR4J(
      P = RainDat$Value,
      Q = FlowDat$Value,
      E = EvapShuffle[,r],
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
  virObsFlowShuffle[,r] <- outputGR4JShuffle$Qsim
}
plotFlowDurationCurve(virObsFlowShuffle, virObsFlow,"withCILimit")
}

dev.off()
