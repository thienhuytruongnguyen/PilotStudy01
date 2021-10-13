#BREE WAS HERE
#Require local WGEN package
install.packages("WGEN",repos = NULL,type = "source")
WD<-getwd()
##Get flow, rain and evap data
RainDat <- read.csv(paste(WD,"/Data/Catchment02/70317_SILO_Rain.csv",sep = ""))
EvapDat <- read.csv(paste(WD,"/Data/Catchment02/70317_SILO_Evap.csv",sep = ""))
FlowDat <- read.csv(paste(WD,"/Data/Catchment02/410730_HRS_FLow.csv",sep = ""))
FlowDat$Value <- (FlowDat$Value*(1E12)/(1.3E14))
####---Add 0.1 to FlowDat
FlowDat_1 <- FlowDat[,2]+0.1

##Fit the WGEN model and Generate some rainfall replicates
RainDatFormat <- WGEN::format_TimeSeries(RainDat)
SimRainList <- WGEN::getSimRain(RainDatFormat, rep = 10, mod = "gama")

##---------------------------------------------------------------------#

library(airGR)
DATA <- data.frame(matrix(NA,nrow = length(RainDat[,1]), ncol = 4))
colnames(DATA) <- c("DatesR","P","Q","E")
DATA[,1] <- RainDat[,1]; DATA[,2] <- RainDat[,2]; DATA[,3] <- FlowDat[,2]; DATA[,4] <- EvapDat[,2]


DATA$DatesR <- strptime(as.character(DATA$DatesR), "%d/%m/%Y")
DATA$DatesR <- format(DATA$DatesR,"%Y-%m-%d")
DATA$DatesR <- as.POSIXlt(DATA$DatesR,tz="",format="%Y-%m-%d")#Format date string

#InputsModel object
InputsModel <- airGR::CreateInputsModel(FUN_MOD = airGR::RunModel_GR4J, DatesR = DATA$DatesR,
                                 Precip = DATA$P, PotEvap = DATA$E)

#RunOptions object
##1.Index Run and WarmUp period
Ind_Run <- seq(which(format(DATA$DatesR, format = "%Y-%m-%d") == "1970-01-01"),
               which(format(DATA$DatesR, format = "%Y-%m-%d") == "2019-02-28"))
Ind_WarmUp <- seq(which(format(DATA$DatesR, format = "%Y-%m-%d") == "1964-01-01"),
               which(format(DATA$DatesR, format = "%Y-%m-%d") == "1969-12-31"))
##2.Run Option
RunOptions <- airGR::CreateRunOptions(FUN_MOD = airGR::RunModel_GR4J,
                               InputsModel = InputsModel, IndPeriod_Run = Ind_Run,
                               IniStates = NULL, IniResLevels = NULL, IndPeriod_WarmUp = Ind_WarmUp)

#InputsCrit object
InputsCrit <- airGR::CreateInputsCrit(FUN_CRIT = airGR::ErrorCrit_NSE, InputsModel = InputsModel, 
                               RunOptions = RunOptions, VarObs = "Q", Obs = DATA$Q[Ind_Run])

#CalibOptions object
CalibOptions <- airGR::CreateCalibOptions(FUN_MOD = airGR::RunModel_GR4J, FUN_CALIB = airGR::Calibration_Michel)

#CALIBRATION
OutputsCalib <- airGR::Calibration_Michel(InputsModel = InputsModel, RunOptions = RunOptions,
                                   InputsCrit = InputsCrit, CalibOptions = CalibOptions,
                                   FUN_MOD = airGR::RunModel_GR4J)
#Get GR4J parameter
Param <- OutputsCalib$ParamFinalR

#Run GR4J
OutputsModel <- airGR::RunModel_GR4J(InputsModel = InputsModel, RunOptions = RunOptions, Param = Param)

#Plot
plot(OutputsModel, Qobs = DATA$Q[Ind_Run])

#------------------------------Global Optimization------------------------------------------------#
lowerGR4J <- rep(-9.99, times = 4)
upperGR4J <- rep(+9.99, times = 4)

#Differential Evolution
optDE <- DEoptim::DEoptim(fn = OptimGR4J,
                          lower = lowerGR4J, upper = upperGR4J,
                          control = DEoptim::DEoptim.control(NP = 40, trace = 10))
#Particle Swarm
optPSO <- hydroPSO::hydroPSO(fn = OptimGR4J,
                             lower = lowerGR4J, upper = upperGR4J,
                             control = list(write2disk = FALSE, verbose = FALSE))

#-------------------------------Multi-objective optimization---------------------------------------#

InputsCrit_inv <- InputsCrit
InputsCrit_inv$transfo <- "inv"
algo <- "caRamel"
optMO <- caRamel::caRamel(nobj = 2,
                          nvar = 4,
                          minmax = rep(TRUE, 2),
                          bounds = matrix(c(lowerGR4J, upperGR4J), ncol = 2),
                          func = MOptimGR4J,
                          popsize = 100,
                          archsize = 100,
                          maxrun = 15000,
                          prec = rep(1.e-3, 2),
                          carallel = FALSE,
                          graph = FALSE)
param_optMO <- apply(optMO$parameters, MARGIN = 1, FUN = function(x) {
  airGR::TransfoParam_GR4J(x, Direction = "TR")
})
RunOptions$Outputs_Sim <- "Qsim"
run_optMO <- apply(optMO$parameters, MARGIN = 1, FUN = function(x) {
  airGR::RunModel_GR4J(InputsModel = InputsModel,
                       RunOptions = RunOptions,
                       Param = x)
}$Qsim)
run_optMO <- data.frame(run_optMO)

x<-BasinObs
ind_t<-seq(which(format(DATA$DatesR, format = "%Y-%m-%d")=="1984-01-01"),
           which(format(DATA$DatesR, format = "%Y-%m-%d")=="2012-12-31"))
x[,2] <- DATA$P[ind_t]
x[,4] <- DATA$E[ind_t]
x[,6] <- DATA$Q[ind_t]
## loading catchment data
data(L0123001)

## preparation of InputsModel object
InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR4J, DatesR = x$DatesR,
                                 Precip = x$P, PotEvap = x$E)

## calibration period selection
Ind_Run <- seq(which(format(x$DatesR, format = "%Y-%m-%d")=="1990-01-01"),
               which(format(x$DatesR, format = "%Y-%m-%d")=="1999-12-31"))

## preparation of RunOptions object
RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR4J, InputsModel = InputsModel,
                               IndPeriod_Run = Ind_Run)

## calibration criterion: preparation of the InputsCrit object
InputsCrit <- CreateInputsCrit(FUN_CRIT = ErrorCrit_NSE, InputsModel = InputsModel,
                               RunOptions = RunOptions, Obs = x$Qmm[Ind_Run])

## preparation of CalibOptions object
CalibOptions <- CreateCalibOptions(FUN_MOD = RunModel_GR4J, FUN_CALIB = Calibration_Michel)

## calibration
OutputsCalib <- Calibration_Michel(InputsModel = InputsModel, RunOptions = RunOptions,
                                   InputsCrit = InputsCrit, CalibOptions = CalibOptions,
                                   FUN_MOD = RunModel_GR4J)

## simulation
Param <- OutputsCalib$ParamFinalR
OutputsModel <- RunModel_GR4J(InputsModel = InputsModel,
                              RunOptions = RunOptions, Param = Param)

## results preview
plot(OutputsModel, Qobs = x$Qmm[Ind_Run])
