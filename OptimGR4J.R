OptimGR4J <- function(ParamOptim) {
  ## Transformation of the parameter set to real space
  RawParamOptim <- airGR::TransfoParam_GR4J(ParamIn = ParamOptim,
                                            Direction = "TR")
  ## Simulation given a parameter set
  OutputsModel <- airGR::RunModel_GR4J(InputsModel = InputsModel,
                                       RunOptions = RunOptions,
                                       Param = RawParamOptim)
  ## Computation of the value of the performance criteria
  OutputsCrit <- airGR::ErrorCrit_RMSE(InputsCrit = InputsCrit,
                                       OutputsModel = OutputsModel,
                                       verbose = FALSE)
  return(OutputsCrit$CritValue)
}