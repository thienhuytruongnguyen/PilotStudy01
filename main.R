#BREE WAS HERE
#Require local WGEN package
install.packages("WGEN",repos = NULL,type = "source")

##Get flow, rain and evap data
RainDat <- read.csv("69041_SILO_Rain.csv")
EvapDat <- read.csv("69041_SILO_Evap.csv")
FlowDat <- read.csv("215004_HRS_Flow.csv")

####---Add 0.1 to FlowDat
FlowDat_1 <- FlowDat[,2]+0.1

##Fit the WGEN model and Generate some rainfall replicates
RainDatFormat <- WGEN::format_TimeSeries(RainDat)
SimRainList <- WGEN::getSimRain(RainDatFormat, rep = 10, mod = "gama")

##---------------------------------------------------------------------#

library(airGR)
DATA <- data.frame(matrix(NA,nrow = length(RainDat[,1]), ncol = 4))
colnames(DATA) <- c("Date_Time","P","Q","E")
DATA[,1] <- RainDat[,1]; DATA[,2] <- RainDat[,2]; DATA[,3] <- FlowDat_1; DATA[,4] <- EvapDat[,2]

DATA$Date_Time <- strptime(as.character(DATA$Date_Time), "%d/%m/%Y")
DATA$Date_Time <- format(DATA$Date_Time,"%Y-%m-%d")
#GR4J model run preparation
Dates <- as.POSIXlt(DATA$Date_Time,tz="",format="%Y-%m-%d")#Format date string

#InputsModel object
InputsModel <- airGR::CreateInputsModel(FUN_MOD = airGR::RunModel_GR4J, DatesR = Dates,
                                 Precip = DATA$P, PotEvap = DATA$E)

#RunOptions object
##1.Index Run and WarmUp period
Ind_Run <- seq(which(format(Dates, format = "%Y-%m-%d") == "1960-01-01"),
               which(format(Dates, format = "%Y-%m-%d") == "2019-02-28"))
Ind_WarmUp <- seq(which(format(Dates, format = "%Y-%m-%d") == "1950-01-01"),
               which(format(Dates, format = "%Y-%m-%d") == "1959-12-31"))
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
