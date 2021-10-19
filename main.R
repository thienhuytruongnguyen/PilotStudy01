##Get work directory
WD<-getwd()

##Get site list
flowSiteList <- read.table("siteList.txt",header=T,sep = "",dec = ".")
weatherSiteList <- getWeatherSiteList()
whichSite <- flowSiteList[3,]##1st row of flow site list

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

##Simulate flow time series
outputGR4J <- runGR4J(paramGR4J)
plot(outputGR4J, Qobs = inputGR4J$Q[paramGR4J[[2]]])

#----------------------------------------------------#

##Fit the WGEN model and Generate some rainfall replicates
RainDatFormat <- format_TimeSeries(RainDat)
rep = 20 #set number of replicates
SimRainList <- getSimRain(RainDatFormat, rep = rep, mod = "gama")
wetday_monthlystats_plot(RainDatFormat,SimRainList)
monthlytotal_stats_plot(RainDatFormat,SimRainList)

##Get Average SimRain
averageSimRain <- getAverageSimRain(SimRainList)
##Fit GR4J model with SimRain
inputGR4J_simRain <- makeInputGR4J(P = averageSimRain, Q = FlowDat$Value, E = EvapDat$Value, Date = FlowDat$Date_Time, A = whichSite$area)
start <- FlowDat[1,1] ; end <- FlowDat[length(FlowDat[,1]),1] #Set start and end of period
warmup <- 24 #months, set warm-up period, airGR prefer more than 12 months
paramGR4J_simRain <- getParamGR4J(inputGR4J = inputGR4J_simRain, start = start, end = end, warmup = warmup)

##Simulate flow time series
outputGR4J_simRain <- runGR4J(paramGR4J_simRain)
plot(outputGR4J_simRain, Qobs = inputGR4J_simRain$Q[paramGR4J[[2]]])

