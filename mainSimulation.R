##--This script is used to run the rainfall-runoff simulation with known GR4J parameters--##
runStart <- Sys.time()

WD<-getwd()

# pdf(
#   file = paste(
#     WD,
#     "/Results&Plots/RainfallStats_FlowDurationCurve_EachSite/",
#     "AnnualMaxvsFDCwithNSE",
#     sep = ""
#   ),
#   width = 8.25,
#   height = 11.75
# )
# 
# par(mfrow=c(2,2))

#Read Site and Silo list info----
flowSiteList <-
  read.table("siteList.txt",
             header = T,
             sep = "",
             dec = ".")
weatherSiteList <- getWeatherSiteList()

catchmentInfo <- read.csv("catchmentInfo.csv")

#Format Station ID for matching data----
weatherSiteList$station <- as.numeric(weatherSiteList$station)
catchmentInfo$NearestSilo <- as.numeric(catchmentInfo$NearestSilo)
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
  
  
  #------------------Rainfall Model----------------------------------
  
  ##Fit the WGEN model and Generate some rainfall replicates
  indRainDate <- makeObsDates(RainDat[,1]) #get index rain day
  RainDatFormat <- format_TimeSeries(RainDat)
  
  rep = 100 #set number of replicates
 
  simRainRep <-
    getSimRain(RainDatFormat,
               rep = rep,
               mod = "gama",
               option = "MoM",
               threshold = 0,
               indRainDate = indRainDate)

  ##Get SimRain Rep
  SimRainList <- makeRainList(simRainRep = simRainRep[[1]], indRainDate = indRainDate)
  
  ##Write rainfall model parameters for the site
  paramWGEN <- data.frame(matrix(NA, ncol = 4, nrow = 12))
  colnames(paramWGEN) <- c("PDW","PWW","a","b")
  
  #Markov Chain Model parameter
  paramWGEN[1:2] <- simRainRep[[2]][1:2]
  
  #Gama distribution parameter
  paramWGEN[3:4] <- simRainRep[[3]]
  
  #Write to .csv file
  write.csv(
    paramWGEN,
    file = paste(
      WD,
      "/WGENmodel_parameter_eachSite/",
      "Site",
      s,
      "_",
      whichSite$station,
      ".csv",
      sep = ""
    ),
    row.names = FALSE
  )
  
  #------------------Get SimFlow Rep----------------------------------

  simFlowRep <-
    getSimFlowRep(simRainRep = simRainRep[[1]], paramGR4J = paramGR4J)
  
  #Get virtual observed flow
  virObsFlow <- outputGR4J$Qsim
  
  
  #Plotting----------------------------------
  plotMaster(WD,
            s,
            whichSite,
            RainDat,
            simRainRep[[1]],
            siloInfo,
            simFlowRep,
            virObsFlow,
            start,
            end,
            paramGR4J,
            SimRainList)
  
}
dev.off()
runEnd <- Sys.time()
runEnd - runStart

obsExceed <- getExceedProb_V2.0(RainDat[paramGR4J[[2]],2])
plot(virObsFDC$flow[14702:19597],obsExceed$flow[14702:19597])
boxplot.ext(virObsFDC$flow)
plot(virObsFDC$Prob,obsExceed$Prob)
par(mfrow=c(3,4))
for (i in 1:12){
  obsExceed <- getExceedProb_V2.0(RainDat[indRainDate$i.mm[[i]],2])
  virObsFDC_Jan <- getExceedProb_V2.0(virObsFlow[indRainDate$i.mm[[i]]])
  plot(obsExceed$flow, virObsFDC_Jan$flow)
}
indRainDate <- makeObsDates(RainDat[paramGR4J[[2]],1])
par(mfrow=c(3,4))
month <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
for (i in 1:12){
  plot(virObsFlow[indRainDate$i.mm[[i]]],RainDat[indRainDate$i.mm[[i]],2],xlab =" ", ylab =" ")
  title(ylab = "Rainfall (mm)", xlab = "Runoff (mm)")
  mtext(paste(month[i]))
}
compareAnnualMaxima(indObsDate = indRainDate, obs = virObsFlow, simRep = simFlowRep)
plotFlowDurationCurve(simFlowRep = simFlowRep, virObsFlow = virObsFlow, option = "withCILimit")
plot(obsExceed,log="y")
plot(virObsFDC, log="y",xlab="", ylab="")
title(ylab = "Flow (mm/day)", xlab="Exceedance Probability")
mtext("Virtual Observed Flow Duration Curve")
