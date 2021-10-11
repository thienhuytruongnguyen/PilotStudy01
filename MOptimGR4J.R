MOptimGR4J <- function(i) {
  if (algo == "caRamel") {
    ParamOptim <- x[i, ]
  }
  ## Transformation of the parameter set to real space
  RawParamOptim <- airGR::TransfoParam_GR4J(ParamIn = ParamOptim,
                                            Direction = "TR")
  ## Simulation given a parameter set
  OutputsModel <- airGR::RunModel_GR4J(InputsModel = InputsModel,
                                       RunOptions = RunOptions,
                                       Param = RawParamOptim)
  ## Computation of the value of the performance criteria
  OutputsCrit1 <- airGR::ErrorCrit_KGE(InputsCrit = InputsCrit,
                                       OutputsModel = OutputsModel,
                                       verbose = FALSE)
  ## Computation of the value of the performance criteria
  OutputsCrit2 <- airGR::ErrorCrit_KGE(InputsCrit = InputsCrit_inv,
                                       OutputsModel = OutputsModel,
                                       verbose = FALSE)
  return(c(OutputsCrit1$CritValue, OutputsCrit2$CritValue))
}