##Get work directory
WD<-getwd()

##Get site list
flowSiteList <- read.table("siteList.txt",header=T,sep = "",dec = ".")
weatherSiteList <- getWeatherSiteList()
x<- flowSiteList[5,]##1st row of flow site list

##Find Nearest Weather station to Site
nearestStation <- getNearestWeatherStation(flowSiteList, weatherSiteList, whichSite = 5)

##Getting rain and PET data from URL
whichStation = list(
  start=x$start, 
  finish=x$finish,
  station=nearestStation,
  format="csv",
  comment="RP",
  username="truonghuythien.nguyen@adelaide.edu.au",
  password="silo"
)
siloData <- getWeatherData(whichStation)

##Get flow, rain and evap data
RainDat <- siloData[,c(2,3)]; colnames(RainDat) = c("Date_Time","Value");
EvapDat <- siloData[,c(2,5)]; colnames(EvapDat) = c("Date_Time","Value")
FlowDat <- read.csv(paste(WD,"/Data/Flow/",x$station,".csv",sep = ""))

##Fit the WGEN model and Generate some rainfall replicates
RainDatFormat <- format_TimeSeries(RainDat)
SimRainList <- getSimRain(RainDatFormat, rep = 10, mod = "gama")
wetday_monthlystats_plot(RainDatFormat,SimRainList)
monthlytotal_stats_plot(RainDatFormat,SimRainList)

##Fit GR4J model
inputGR4J <- makeInputGR4J(P = RainDat$Value, Q = FlowDat$Value, E = EvapDat$Value, Date = FlowDat$Date_Time, A = x$area)
paramGR4J <- getParamGR4J(inputGR4J = inputGR4J)

##Simulate flow time series
outputGR4J <- runGR4J(paramGR4J)
plot(outputGR4J, Qobs = inputGR4J$Q[paramGR4J[[2]]])




