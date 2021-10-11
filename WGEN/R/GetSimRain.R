#' Get simulated rainfall
#'
#' @param occur.param MC model parameters
#' @param amount.param Amount model parameters
#' @param rep Number of desired replicates
#' @param mod Specified "gama" or "expo" distribution
#'
#' @return TS of simulated rainfall
#'
#' @export
getSimRain <- function(obs.data, rep = 100, mod = "gama"){
  rain.data <- format_TimeSeries(obs.data) #formatting the obs data

  ##Declaring model parameters objects
  occur.param <-data.frame()
  amount.param <-data.frame()

  ##Calibrating model parameters
  occur.param <- fitMCModel(rain.data)
  amount.param <- fitAmountModel(rain.data,mod)

  ##Simulating
  ls.month.sim <-Amount_model(occur.param, amount.param, rep, obs.data)

  return(ls.month.sim)
}
