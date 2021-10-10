library(hydromad)library(airGR)

FlowData <- read.csv("streamflow.csv")
RainData <- read.csv("Rainfall.csv")
PetData <- read.csv("PET.csv")

DATA <- data.frame(matrix(NA,nrow = length(RainData[,1]), ncol = 4))
colnames(DATA) <- c("Date_time","P","Q","E")
DATA[,1] <- RainData[,1]; DATA[,2] <- RainData[,2]; DATA[,3] <- FlowData[,2]; DATA[,4] <- PetData[,2]



#GR4J model run preparation
Dates <- as.POSIXlt(DATA$Date_time,tz="",format="%Y-%m-%d")#Format date string

#InputsModel object
InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR4J, DatesR = Dates,
                                 Precip = DATA$P, PotEvap = DATA$E)

#RunOptions object
 ##1.Index Run
Ind_Run <- seq(which(format(Dates, format = "%Y-%m-%d") == "1994-01-01"), #start from row 550 of the input ts
               which(format(Dates, format = "%Y-%m-%d") == "2019-02-28"))

 ##2.Run Option
RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR4J,
                               InputsModel = InputsModel, IndPeriod_Run = Ind_Run,
                               IniStates = NULL, IniResLevels = NULL, IndPeriod_WarmUp = NULL)

#InputsCrit object
InputsCrit <- CreateInputsCrit(FUN_CRIT = ErrorCrit_NSE, InputsModel = InputsModel, 
                               RunOptions = RunOptions, VarObs = "Q", Obs = DATA$Q[550:9739])#!!Same length with run index

#CalibOptions object
CalibOptions <- CreateCalibOptions(FUN_MOD = RunModel_GR4J, FUN_CALIB = Calibration_Michel)

#CALIBRATION
OutputsCalib <- Calibration_Michel(InputsModel = InputsModel, RunOptions = RunOptions,
                                   InputsCrit = InputsCrit, CalibOptions = CalibOptions,
                                   FUN_MOD = RunModel_GR4J)
