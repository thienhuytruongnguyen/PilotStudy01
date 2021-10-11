#' Calibrate the 2-state Markov Cahin model parameters PWW and PDW
#'
#' @param obs.data observed daily TS of rainfall
#'
#' @return MC model parameters PWW and PDW
#'
#' @export
fitMCModel<- function
(obs.data ##<<Formatted observed data
){
  ##Declare a dataframe to store the calibrated parameters
  param <- data.frame(matrix(nrow = 12,ncol = 4))

  for (i in 1:12){##Loop for each month

    x <- occurr_stats(obs.data[obs.data$month==i,2]) #x is a temporary dataframe to store the outcome of function occurr_stats
    param[i,] <- x[1,] #Store the occurrence parameters set for each month

  }

  return(param)

}

