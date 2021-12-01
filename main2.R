##--This script is used to run the rainfall-runoff simulation with known GR4J parameters--##
runStart <- Sys.time()

WD<-getwd()

pdf(
  file = paste(
    WD,
    "/Results&Plots/RainfallStats_FlowDurationCurve_EachSite/",
    "AnnualMaxvsFDCwithNSE",
    sep = ""
  ),
  width = 8.25,
  height = 11.75
)

par(mfrow=c(2,2))

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
  rep = 100 #set number of replicates
  SimRainList <-
    getSimRain(RainDat,
               rep = rep,
               mod = "gama",
               option = "MoM",
               threshold = 0)
  
  ##Get SimRain Rep
  simRainRep <- getSimRainRep(SimRainList[[1]])
  
  ##Write rainfall model parameters for the site
  paramWGEN <- data.frame(matrix(NA, ncol = 4, nrow = 12))
  colnames(paramWGEN) <- c("PDW","PWW","a","b")
  
  #Markov Chain Model parameter
  paramWGEN[1:2] <- SimRainList[[2]][1:2]
  
  #Gama distribution barameter
  paramWGEN[3:4] <- SimRainList[[3]]
  
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
    getSimFlowRep(simRainRep = simRainRep, paramGR4J = paramGR4J)
  
  #Get virtual observed flow
  virObsFlow <- outputGR4J$Qsim
  
  #--------------------Calculate Efficiency--------------------
  
  ##get annual maxima of rainfall-----------------------
  indRainDate <- makeObsDates(RainDat$Date_Time[paramGR4J[[2]]]) ##Get date index
  
  #obs annual maxima
  obsAnnualMax <-
    getAnnualMaxima(indObsDate = indRainDate, value = RainDat$Value[paramGR4J[[2]]])
  
  #sim annual maxima
  simAnnualMaxima <-
    data.frame(matrix(
      NA,
      nrow = length(obsAnnualMax),
      ncol = ncol(simRainRep)
    ))
  
  for (i in 1: ncol(simAnnualMaxima)){
    simAnnualMaxima[, i] <-
      getAnnualMaxima(indObsDate = indRainDate, value = simRainRep[paramGR4J[[2]], i])
  }
  
  # ##get NSE annual maxima----
  # NSE_annualMaxima <- getNSE(obs = obsAnnualMax, sim = simAnnualMaxima)
  # NSE_annualMaxima <- format(NSE_annualMaxima,digits=3)
  # ##get NSE flow duration curve----
  # NSE_flowDurationCurve <- getNSE_FDC(obs = virObsFlow, sim = simFlowRep)
  # NSE_flowDurationCurve <- format(NSE_flowDurationCurve,digits=3)
  
  #Plotting----------------------------------
  pdf(
    file = paste(
      WD,
      "/Results&Plots/RainfallStats_FlowDurationCurve_EachSite/",
      "Site",
      s,
      "_",
      whichSite$station,
      ".pdf",
      sep = ""
    ),
    width = 8.25,
    height = 11.75
  )
  
  par(mfrow=c(2,2))
  
  ##Plot return interval for annual maximum rainfall----
  
  #call plot function
  compareAnnualMaxima(indObsDate = indRainDate, obs = RainDat$Value, simRep = simRainRep)
  
  #get station info
  siloStation <-
    paste(
      " Daily Annual Maxima (Rainfall)",
      "\n",
      " ID:",
      siloInfo$station,
      "Name:",
      as.character(siloInfo$name),
      "\n",
      "Juradiction:",
      siloInfo$state,
      "Lat:",
      siloInfo$latitude,
      "Long:",
      siloInfo$longitude,
      "Elevation:",
      siloInfo$elevation,
      sep = " "
    )
  
  mtext(siloStation, 3, 0, cex=0.5, adj=0, padj = -0.3)
  
  ##Plot Flow Duration Curve for flow depth----
  #Call plot function
  plotFlowDurationCurve(simFlowRep = simFlowRep, virObsFlow = virObsFlow, option = "withCILimit")
  #get outlet info
  outlet <-
    paste(
      " ID:",
      whichSite[1],
      "Juradiction:",
      whichSite[4],
      "Area:",
      whichSite[5],
      "km^2",
      "\n",
      "Lat:",
      whichSite[2],
      "Long:",
      whichSite[3],
      "From:",
      start,
      "To:",
      end,
      sep = " "
    )
  
  mtext(outlet, 3, 0, cex=.5, adj=0, padj=-.3)
  
   ##Plot return interval for annual maximumflow----
  
   indFlowDate <- makeObsDates(RainDat$Date_Time[paramGR4J[[2]]]) #Get date index'
  
   #call plot function
   compareAnnualMaxima(indObsDate = indFlowDate, obs = virObsFlow, simRep = simFlowRep)
   
   mtext("Daily annual maxima (Flow)")

   ##Plot simulated vs observed rainfall stats----
   
   RainDatFormat <- format_TimeSeries(RainDat)
   
   #wet day amount stats
   wetday_monthlystats_plot(RainDatFormat, SimRainList[[1]], type = "boxplot")
   
   #Monthly total stats
   monthlytotal_stats_plot(RainDatFormat, SimRainList[[1]], type = "boxplot")
   
   #mean 3-5 day total monthly
   compareMean3dayTotal(obs = RainDat$Value, sim = simRainRep, indObsDate = indRainDate)
   compareMean5dayTotal(obs = RainDat$Value, sim = simRainRep, indObsDate = indRainDate)
   
   #Dry Spell Plot
   month <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
   
   maxSpell <- 10
   
   for (m in 1:12){
     compareDrySpell(
       monthlyObsRain = RainDat$Value[indRainDate$i.mm[[m]]],
       monthlySimRain = simRainRep[indRainDate$i.mm[[m]], ],
       maxSpell = maxSpell
     )
     mtext(paste("Dry Spell: ", paste(month[m]), sep=""))
     
     #Wet Spell Plot
     compareWetSpell(
       monthlyObsRain = RainDat$Value[indRainDate$i.mm[[m]]],
       monthlySimRain = simRainRep[indRainDate$i.mm[[m]], ],
       maxSpell = maxSpell
     )
     mtext(paste("Wet Spell: ", paste(month[m]), sep=""))
   }
   
  
  dev.off()
}


dev.off()
runEnd <- Sys.time()
runEnd - runStart

#Random shuffle----
simRainListShuffle <- SimRainList
for (m in 1:12){
  for (r in 1:100){
    simRainListShuffle[[1]][[m]][r,] <- randomShuffle(simRainListShuffle[[1]][[m]][r,],10)
  }
}

simRainRepShuffle <- getSimRainRep(simRainListShuffle[[1]])

simFlowRepShuffle <-
  getSimFlowRep(simRainRep = simRainRepShuffle, paramGR4J = paramGR4J)

plotFlowDurationCurve(simFlowRep = simFlowRep, virObsFlow = virObsFlow, option = "withCILimit")
plotFlowDurationCurve(simFlowRep = simFlowRepShuffle, virObsFlow = virObsFlow, option = "withCILimit")
compareAnnualMaxima(indObsDate = indFlowDate, obs = virObsFlow, simRep = simFlowRep)
compareAnnualMaxima(indObsDate = indFlowDate, obs = virObsFlow, simRep = simFlowRepShuffle)
#wet day amount stats
wetday_monthlystats_plot(RainDatFormat, simRainListShuffle[[1]], type = "boxplot")

#Monthly total stats
monthlytotal_stats_plot(RainDatFormat, SimRainList[[1]], type = "boxplot")

#mean 3-5 day total monthly
compareMean3dayTotal(obs = RainDat$Value, sim = simRainRep, indObsDate = indRainDate)
compareMean5dayTotal(obs = RainDat$Value, sim = simRainRep, indObsDate = indRainDate)


#----Using SCE Optimiser----

iniThetaDataframe <-
  read.csv(file = paste(WD, "/WGENmodel_parameter_eachSite/Site1_410730.csv", sep = ""))
iniTheta <- rep(0,48)
iniTheta[1:12]<-iniThetaDataframe[1:12,1]; iniTheta[13:24]<-iniThetaDataframe[1:12,2]
iniTheta[25:36]<-iniThetaDataframe[1:12,3]; iniTheta[37:48]<-iniThetaDataframe[1:12,4]
iniTheta <- theta_Trial02
lowerTheta <- rep(0,48)
lowerTheta[1:24] = 0; lowerTheta[25:48] = 1e-10
upperTheta <- rep(0,48)
upperTheta[1:24] = 1; upperTheta[25:48] = Inf

opt <- hydromad::SCEoptim(FUN = SSE_FlowDurationCurve,
                   par = iniTheta,
                   obsRain = RainDat,
                   paramGR4J = paramGR4J,
                   virObsFlow = virObsFlow,
                   lower = lowerTheta,
                   upper = upperTheta)
Sys.time()



theta_Trial04 <- opt$par

##Collecting results from trial
trial <- list(iniTheta, theta_Trial01, theta_Trial02, theta_Trial03,theta_Trial04,theta_Trial05)

for (t in 1:6){
  theta <- trial[[t]]
  
  simFlowRep_Opt <- getSimFlowRep_Opt(theta = theta,
                                      paramGR4J = paramGR4J,
                                      rep = 100,
                                      obsRain = RainDat)
  
  
  simRainRep_Opt <- getSimRainRep(simFlowRep_Opt[[2]])
  indRainDate <- makeObsDates(RainDat$Date_Time[paramGR4J[[2]]]) ##Get date index
  
  pdf(
    file = paste(
      WD,
      "/Results&Plots/SCEoptim/",
      "Site",
      s,
      "_",
      whichSite$station,
      "_",
      "Trial",
      t,
      ".pdf",
      sep = ""
    ),
    width = 8.25,
    height = 11.75
  )
  
  par(mfrow=c(2,2))
  
  compareAnnualMaxima(indObsDate = indRainDate, obs = RainDat$Value, simRep = simRainRep_Opt)
  
  siloStation <-
    paste(
      " Daily Annual Maxima (Rainfall)",
      "\n",
      " ID:",
      siloInfo$station,
      "Name:",
      as.character(siloInfo$name),
      "\n",
      "Juradiction:",
      siloInfo$state,
      "Lat:",
      siloInfo$latitude,
      "Long:",
      siloInfo$longitude,
      "Elevation:",
      siloInfo$elevation,
      sep = " "
    )
  
  mtext(siloStation,
        3,
        0,
        cex = 0.5,
        adj = 0,
        padj = -0.3)
  
  plotFlowDurationCurve(simFlowRep = simFlowRep_Opt[[1]],
                        virObsFlow = virObsFlow,
                        option = "withCILimit")
  
  outlet <-
    paste(
      " Flow Duration Curve",
      "\n",
      " ID:",
      whichSite[1],
      "Juradiction:",
      whichSite[4],
      "Area:",
      whichSite[5],
      "km^2",
      "\n",
      "Lat:",
      whichSite[2],
      "Long:",
      whichSite[3],
      "From:",
      start,
      "To:",
      end,
      sep = " "
    )
  
  mtext(outlet,
        3,
        0,
        cex = .5,
        adj = 0,
        padj = -.3)
  
  compareAnnualMaxima(indObsDate = indFlowDate,
                      obs = virObsFlow,
                      simRep = simFlowRep_Opt[[1]])
  
  mtext("Daily annual maxima (Flow)")
  
  wetday_monthlystats_plot(RainDatFormat, simFlowRep_Opt[[2]], type = "boxplot")
  monthlytotal_stats_plot(RainDatFormat, simFlowRep_Opt[[2]], type = "boxplot")
  compareMean3dayTotal(obs = RainDat$Value, sim = simFlowRep_Opt[[3]], indObsDate = indRainDate)
  compareMean5dayTotal(obs = RainDat$Value, sim = simFlowRep_Opt[[3]], indObsDate = indRainDate)
  #Dry Spell Plot
  month <-
    c("Jan",
      "Feb",
      "Mar",
      "Apr",
      "May",
      "Jun",
      "Jul",
      "Aug",
      "Sep",
      "Oct",
      "Nov",
      "Dec")
  
  maxSpell <- 10
  
  for (m in 1:12){
    compareDrySpell(
      monthlyObsRain = RainDat$Value[indRainDate$i.mm[[m]]],
      monthlySimRain = simRainRep_Opt[indRainDate$i.mm[[m]], ],
      maxSpell = maxSpell
    )
    mtext(paste("Dry Spell: ", paste(month[m]), sep=""))
    
    #Wet Spell Plot
    compareWetSpell(
      monthlyObsRain = RainDat$Value[indRainDate$i.mm[[m]]],
      monthlySimRain = simRainRep_Opt[indRainDate$i.mm[[m]], ],
      maxSpell = maxSpell
    )
    mtext(paste("Wet Spell: ", paste(month[m]), sep=""))
  }
  
  dev.off()
}



  paramMC <- data.frame(matrix(NA,12,2))
  paramMC[,1] <- iniTheta[1:12]; paramMC[,2] <- iniTheta[13:24]
  paramAmount <- data.frame(matrix(NA,12,2))
  paramAmount[,1] <- iniTheta[25:36]; paramAmount[,2] <- iniTheta[37:48]

write.csv(paramMC,"occur.csv")
write.csv(paramAmount,"amount.csv")
opt <- optim(fn = SSE_FlowDurationCurve,
                          par = iniTheta,
                          obsRain = RainDat,
                          paramGR4J = paramGR4J,
                          virObsFlow = virObsFlow,
                          lower = lowerTheta,
                          upper = upperTheta, method = "L-BFGS-B")



