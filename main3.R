##--This script is used to run the rainfall-runoff simulation with known GR4J parameters--
## add option for Method of moment for amount model
## Status: MoM does not work
runStart <- Sys.time()
WD<-getwd()
pdf(file = paste(WD,"/FDC_results/Fixed/","AnnualMaxvsFDCwithNSE",sep = ""), 
    width = 8.25, height = 11.75)
par(mfrow=c(2,2))

#Read Site and Silo list info
flowSiteList <- read.table("siteList.txt",header=T,sep = "",dec = ".")
weatherSiteList <- getWeatherSiteList()
catchmentInfo <- read.csv("catchmentInfo.csv")

#Format Station ID for matching data
weatherSiteList$station <- as.numeric(weatherSiteList$station)
catchmentInfo$NearestSilo <- as.numeric(catchmentInfo$NearestSilo)


for (i in 1:25) {
  
  # pdf(file = paste(WD,"/FDC_results/Fixed/","Site",i,".pdf",sep = ""), 
  #     width = 8.25, height = 11.75)
  # par(mfrow=c(2,2))
  
  #------------------Get Data----------------------------------#
  
  ##Get site list
  whichSite <- flowSiteList[i,]##1st row of flow site list
  
  ##Get Nearest Weather station to Site
  siloInfo <- weatherSiteList[which(weatherSiteList$station==catchmentInfo$NearestSilo[i]),]
  
  ##Getting rain and PET data from URL
  siloData <- getWeatherData(nearestStation = siloInfo$station, start = whichSite$start, finish = whichSite$finish)
  
  ##Get flow, rain and evap data
  RainDat <- siloData[,c(2,3)]; colnames(RainDat) = c("Date_Time","Value");
  EvapDat <- siloData[,c(2,5)]; colnames(EvapDat) = c("Date_Time","Value")
  FlowDat <- read.csv(paste(WD,"/Data/Flow/",whichSite$station,".csv",sep = ""))
  
  #------------------Rainfall Model MLE----------------------------------#
  
  ##Fit the WGEN model and Generate some rainfall replicates
  rep = 100 #set number of replicates
  SimRainList_MLE <- getSimRain(RainDat, rep = rep, mod = "gama", option = "MLE")
  
  ##Get SimRain Rep
  simRainRep_MLE<- getSimRainRep(SimRainList_MLE)
  
  #------------------Rainfall Model MoM----------------------------------#
  ##Fit the WGEN model and Generate some rainfall replicates
  rep = 100 #set number of replicates
  SimRainList_MoM <- getSimRain(RainDat, rep = rep, mod = "gama", option = "MoM")
  
  ##Get SimRain Rep
  simRainRep_MoM<- getSimRainRep(SimRainList_MoM)
  
  #------------------Runoff Model----------------------------------#

  ##Set GR4J model parameter
  inputGR4J <- makeInputGR4J(P = RainDat$Value, Q = FlowDat$Value, E = EvapDat$Value, Date = FlowDat$Date_Time, A = whichSite$area)
  start <- FlowDat[1,1] ; end <- FlowDat[length(FlowDat[,1]),1] #Set start and end of period
  warmup <- 24 #months, set warm-up period, airGR prefer more than 12 months
  paramGR4J <- getParamGR4J(inputGR4J = inputGR4J, start = start, end = end, warmup = warmup, parameter = "known", knownParam = catchmentInfo[i,c(9,10,11,12)])

  ##Run GR4J model
  outputGR4J <- runGR4J(paramGR4J)
  effiReport <- getEffiReport(outputGR4J = outputGR4J, inputGR4J = inputGR4J, start = start, end = end, warmup = warmup)

  ##Get SimFlow Rep
  simFlowRep <- getSimFlowRep(simRainRep = simRainRep_MoM, paramGR4J = paramGR4J)

  #Get virtual observed flow
  virObsFlow <- outputGR4J$Qsim

  #--------------------calculate Efficiency--------------------#

  ##get annual maxima of rainfall-----------------------
  indRainDate <- makeObsDates(RainDat$Date_Time) ##Get date index
  #obs annual maxima
  obsAnnualMax <- getAnnualMaxima(indObsDate = indRainDate, value = RainDat$Value)
  #sim annual maxima
  simAnnualMaxima <- data.frame(matrix(NA, nrow = length(obsAnnualMax), ncol = ncol(simRainRep)))
  for (i in 1: ncol(simAnnualMaxima)){
    simAnnualMaxima[,i] <- getAnnualMaxima(indObsDate = indRainDate, value = simRainRep[,i])
  }
  #get NSE annual maxima
  NSE_annualMaxima <- getNSE(obs = obsAnnualMax, sim = simAnnualMaxima)
  NSE_annualMaxima <- format(NSE_annualMaxima,digits=3)
  ##get NSE flow duration curve----------------------
  NSE_flowDurationCurve <- getNSE_FDC(obs = virObsFlow, sim = simFlowRep)
  NSE_flowDurationCurve <- format(NSE_flowDurationCurve,digits=3)

  #------------------Plotting----------------------------------#

  ##Plot return interval for annual maximum rainfall

  #call plot function
  compareAnnualMaxima(indObsDate = indRainDate, obs = RainDat$Value, simRep = simRainRep_MoM)

  siloStation <- paste(" ID:",siloInfo$station, "Name:", as.character(siloInfo$name), "\n", "Juradiction:", siloInfo$state, "Lat:", siloInfo$latitude, "Long:", siloInfo$longitude, "Elevation:", siloInfo$elevation, "NSE:",NSE_annualMaxima, sep = " ")
  mtext(siloStation, 3, 0, cex=0.5, adj=0, padj = -0.3)

  ##Plot Flow Duration Curve for flow depth
  #Call plot function
  plotFlowDurationCurve(simFlowRep = simFlowRep, virObsFlow = virObsFlow, option = "withCILimit")

  outlet <- paste(" ID:",whichSite[1],"Juradiction:",whichSite[4],"Area:",
                  whichSite[5],"km^2","\n","Lat:",whichSite[2],"Long:",whichSite[3],"From:",start,"To:",end,"NSE:",NSE_flowDurationCurve,sep = " ")
  mtext(outlet, 3, 0, cex=.5, adj=0, padj=-.3)

   ##Plot return interval for annual maximum flow
   indFlowDate <- makeObsDates(RainDat$Date_Time[paramGR4J[[2]]]) #Get date index'
  #call plot function
   compareAnnualMaxima(indObsDate = indFlowDate, obs = virObsFlow, simRep = simFlowRep)

   #Plot simulated vs observed rainfall stats
   RainDatFormat <- format_TimeSeries(RainDat)
  #MLE-----------
   wetday_monthlystats_plot(RainDatFormat,SimRainList_MLE, type = "boxplot") #plot monthly wet day stats, type of plot can be "boxplot" or "errorbar"
   monthlytotal_stats_plot(RainDatFormat,SimRainList_MLE, type = "boxplot") #plot monthly total stats, type of plot can be "boxplot" or "errorbar"
   #MoM-----------
   wetday_monthlystats_plot(RainDatFormat,SimRainList_MoM, type = "boxplot") #plot monthly wet day stats, type of plot can be "boxplot" or "errorbar"
   monthlytotal_stats_plot(RainDatFormat,SimRainList_MoM, type = "boxplot") #plot monthly total stats, type of plot can be "boxplot" or "errorbar"
  #dev.off()
}
#########
dev.off()
runEnd <- Sys.time()
runEnd - runStart


