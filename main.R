#Require local WGEN package
install.packages("WGEN",repos = NULL,type = "source")
library(WGEN)

##Get flow, rain and evap data
RainDat <- read.csv("69049_SILO_Rain.csv")
EvapDat <- read.csv("69049_SILO_Rain.csv")
FlowDat <- read.csv("215004_HRS_Flow.csv")

##Fit the WGEN model and Generate some rainfall replicates

library(WGEN)
RainDatFormat <- format_TimeSeries(RainDat)
SimRainList <- getSimRain(RainDatFormat, rep = 10, mod = "expo")

