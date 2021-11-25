##-----This script is used to fit the WGEN model and the GR4J model and run the rainfall-runoff simulation----##

##Get work directory
runStart <- Sys.time()
WD<-getwd()
pdf(file = paste(WD,"/FDC_results/Fixed/","AnnualMaxvsFDC.pdf",sep = ""), 
    width = 8.25, height = 11.75)
par(mfrow=c(2,2))

catchmentParam <- data.frame(matrix(NA,nrow = 25, ncol = 13))
colnames(catchmentParam) <- c("ID","Lat","Long","Juradiction","Area(Km^2)","Start","End","x1","x2","x3","x4","NSE","NearestSilo")
flowSiteList <- read.table("siteList.txt",header=T,sep = "",dec = ".")
weatherSiteList <- getWeatherSiteList()
for (i in 1:25) {
  ##Get site list
  whichSite <- flowSiteList[i,]##1st row of flow site list
  catchmentParam[i,1:7] <- whichSite[,1:7]
  ##Find Nearest Weather station to Site
  nearestStation <- getNearestWeatherStation(flowSiteList, weatherSiteList, whichSite = whichSite$ID)
  catchmentParam[i,13] <- nearestStation[[1]]
  ##Getting rain and PET data from URL
  siloData <- getWeatherData(nearestStation = nearestStation[[1]], start = whichSite$start, finish = whichSite$finish)
  
  ##Get flow, rain and evap data
  RainDat <- siloData[,c(2,3)]; colnames(RainDat) = c("Date_Time","Value");
  EvapDat <- siloData[,c(2,5)]; colnames(EvapDat) = c("Date_Time","Value")
  FlowDat <- read.csv(paste(WD,"/Data/Flow/",whichSite$station,".csv",sep = ""))
  
  ##Fit GR4J model
  inputGR4J <- makeInputGR4J(P = RainDat$Value, Q = FlowDat$Value, E = EvapDat$Value, Date = FlowDat$Date_Time, A = whichSite$area)
  start <- FlowDat[1,1] ; end <- FlowDat[length(FlowDat[,1]),1] #Set start and end of period
  warmup <- 24 #months, set warm-up period, airGR prefer more than 12 months
  paramGR4J <- getParamGR4J(inputGR4J = inputGR4J, start = start, end = end, warmup = warmup)
  
  ##Get NSE report and store parameters
  catchmentParam[i,8:11] <- format(round(paramGR4J[[1]],3)) 
  
  ##Run GR4J model
  outputGR4J <- runGR4J(paramGR4J)
  effiReport <- getEffiReport(outputGR4J = outputGR4J, inputGR4J = inputGR4J, start = start, end = end, warmup = warmup)
  catchmentParam[i,12] <- format(round(effiReport$CritValue,3))
  
  #----------------------------------------------------#
  
  ##Fit the WGEN model and Generate some rainfall replicates
  rep = 100 #set number of replicates
  SimRainList <- getSimRain(RainDat, rep = rep, mod = "gama")
  
  ##Get SimRain Rep
  simRainRep<- getSimRainRep(SimRainList)
  
  ##Get SimFlow Rep
  simFlowRep <- getSimFlowRep(simRainRep = simRainRep, paramGR4J = paramGR4J)
  
  #Get virtual observed flow
  virObsFlow <- outputGR4J$Qsim
  
  #Plot return interval
  indObsDate <- makeObsDates(RainDat$Date_Time)
  station <- paste(" ID:",nearestStation[[2]]$station, "Name:", as.character(nearestStation[[2]]$name), "\n", "Juradiction:", nearestStation[[2]]$state, "Lat:", nearestStation[[2]]$latitude, "Long:", nearestStation[[2]]$longitude, "Elevation:", nearestStation[[2]]$elevation)
  
  compareAnnualMaxima(indObsDate = indObsDate, obsRain = RainDat$Value, simRainRep = simRainRep)
  
  mtext(station, 3, 0, cex=0.5, adj=0, padj = -0.3)
  
  #Plot FDC
  outlet <- paste(" ID:",whichSite[1],"Juradiction:",whichSite[4],"Area:",
                   whichSite[5],"km^2","\n","Lat:",whichSite[2],"Long:",whichSite[3],"From:",start,"To:",end,sep = " ")
  
  plotFlowDurationCurve(simFlowRep = simFlowRep, virObsFlow = virObsFlow, option = "withCILimit")
  
  mtext(outlet, 3, 0, cex=.5, adj=0, padj=-.3)
  
}
##Write parameters table with site info
# ttable <- gridExtra::tableGrob(catchmentParam, rows=NULL, theme = gridExtra::ttheme_default(base_size = 8.5))
# grid::grid.draw(ttable)
write.csv(catchmentParam, file = "catchmentInfo.csv")
    
#########
dev.off()
runEnd <- Sys.time()
runEnd - runStart
#########
    


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

plot(outputGR4J_simRain, Qobs = inputGR4J$Q[paramGR4J_simRain[[2]]])+#plot comparison obs vs simflow_simRain
mtext("Obseved Flow vs Simulated Flow with Simulated Rain", side = 4, line = 0.5)

#Plot FDC
plotFlowDurationCurve(simFlowRep = simFlowRep, virObsFlow = virObsFlow, option = "RepVSVirObs")

##########################################
WD<-getwd()

##Get site list
flowSiteList <- read.table("siteList.txt",header=T,sep = "",dec = ".")
weatherSiteList <- getWeatherSiteList()
whichSite <- flowSiteList[1,]##1st row of flow site list

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

#Get virtual observed flow
virObsFlow <- outputGR4J$Qsim

#Get WGEN parameters
paramMC <- fitMCModel(RainDat)
paramAmount <- fitAmountModel(RainDat,mod = "gama")
colnames(paramAmount) <- c("alpha","beta")

#make WGEM parameter
paramMC
paramAmount$alpha[c(6,7,8)] <- paramAmount$alpha[c(6,7,8)]*8
paramAmount$alpha[c(12,1,2)] <- paramAmount$alpha[c(12,1,2)]/8
paramAmount$alpha[c(3,4,5)] <- paramAmount$alpha[c(3,4,5)]*8
paramAmount$alpha[c(9,10,11)] <- paramAmount$alpha[c(9,10,11)]/8
paramAmount$beta <- paramAmount$beta

#plot FDC with changing WGEN parameters
simRaintoFDCModel(paramMC = paramMC,
                  paramAmount = paramAmount,
                  obs.data = RainDat,
                  rep = 10,
                  paramGR4J = paramGR4J,
                  virObsFlow = virObsFlow)



