##Get work directory
WD<-getwd()

##Get site list
flowSiteList <- read.table("siteList.txt",header=T,sep = "",dec = ".")
weatherSiteList <- getWeatherSiteList()
whichSite <- flowSiteList[24,]##1st row of flow site list

##Find Nearest Weather station to Site
nearestStation <- getNearestWeatherStation(flowSiteList, weatherSiteList, whichSite = whichSite$ID)

##Getting rain and PET data from URL
siloData <- getWeatherData(whichStation, start = whichSite$start, finish = whichSite$finish)

##Get flow, rain and evap data
RainDat <- siloData[,c(2,3)]; colnames(RainDat) = c("Date_Time","Value");
EvapDat <- siloData[,c(2,5)]; colnames(EvapDat) = c("Date_Time","Value")
FlowDat <- read.csv(paste(WD,"/Data/Flow/",whichSite$station,".csv",sep = ""))
FlowDat[,2] <- FlowDat[,2] + 0.01
##Fit GR4J model
inputGR4J <- makeInputGR4J(P = RainDat$Value, Q = FlowDat$Value, E = EvapDat$Value, Date = FlowDat$Date_Time, A = whichSite$area)
start <- FlowDat[1,1] ; end <- FlowDat[length(FlowDat[,1]),1] #Set start and end of period
warmup <- 24 #months, set warm-up period, airGR prefer more than 12 months
paramGR4J <- getParamGR4J(inputGR4J = inputGR4J, start = start, end = end, warmup = warmup)

##Simulate flow time series
outputGR4J <- runGR4J(paramGR4J)

obsEffiReport <- getEffiReport(outputGR4J = outputGR4J, inputGR4J = inputGR4J, start = start, end = end, warmup = warmup) #Efficiency report
#----------------------------------------------------#

##Fit the WGEN model and Generate some rainfall replicates
RainDatFormat <- format_TimeSeries(RainDat)
rep = 20 #set number of replicates
SimRainList <- getSimRain(RainDatFormat, rep = rep, mod = "gama")

##Get Average SimRain
SimRain <- getSimRainRep(SimRainList)

##make input and param list with average simRain
inputGR4J_simRain <- makeInputGR4J(P = SimRain[,3], Q = FlowDat$Value, E = EvapDat$Value, Date = FlowDat$Date_Time, A = whichSite$area)

paramGR4J_simRain <- paramGR4J
paramGR4J_simRain[[3]]$Precip <- SimRain[,3]

##Simulate flow time series with average simRain
outputGR4J_simRain <- runGR4J(paramGR4J_simRain)

simEffiReport <- getEffiReport(outputGR4J = outputGR4J_simRain, inputGR4J = inputGR4J_simRain, start = start, end = end, warmup = warmup) #Efficiency report

##Plot

#Plot simulated vs observed rainfall stats
wetday_monthlystats_plot(RainDatFormat,SimRainList, type = "boxplot") #plot monthly wet day stats, type of plot can be "boxplot" or "errorbar"
monthlytotal_stats_plot(RainDatFormat,SimRainList, type = "boxplot") #plot monthly total stats, type of plot can be "boxplot" or "errorbar"

#plot simulated vs observed flow stats
plot(outputGR4J, Qobs = inputGR4J$Q[paramGR4J[[2]]])+#plot comparison obs vs simflow_obsRain
mtext("Obseved Flow vs Simulated Flow with Observed Rain", side = 4, line = 0.5)

plot(outputGR4J_simRain, Qobs = inputGR4J_simRain$Q[paramGR4J_simRain[[2]]])+#plot comparison obs vs simflow_simRain
mtext("Obseved Flow vs Simulated Flow with Simulated Rain", side = 4, line = 0.5)

virtualObsExceedProb <- getExceedProb(outputGR4J$Qsim)
plot(simExceedProb,log="y",type="l", col = rgb(red = 1, green = 0, blue = 0, alpha = 0.2),ylim=c(min(simExceedProb$Flow),max(simExceedProb)),lwd=2)
lines(virtualObsExceedProb,col="blue",lwd=4)
simExceedProb <- getExceedProb(outputGR4J_simRain$Qsim)
