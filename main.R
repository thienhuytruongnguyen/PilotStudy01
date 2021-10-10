##Get flow, rain and evap data
RainDat <- read.csv("69049_SILO_Rain.csv")
EvapDat <- read.csv("69049_SILO_Rain.csv")
FlowDat <- read.csv("215004_HRS_Flow.csv")

##Fit the WGEN model and Generate some rainfall replicates
#>install.packages("C:/Research/Programming/PilotStudy01/Packages/TSformat",repos = NULL, type = "source")
library(FormatTimeSeries)

##Update from T laptop
##Delete T laptop line