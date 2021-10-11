#' The Amount model
#'
#' @param occur.param MC model parameters
#' @param amount.param Amount model parameters
#' @param rep Number of desired replicates
#' @param obs.data Formatted obs data
#'
#' @return TS of simulated rainfall
#'
#' @export
Amount_model <- function
(occur.param,##<<Markov model parameters, PWW and PDW
 amount.param,##<<Amount model parameters, \lambda for exponential distrib or \alpha and \beta for gamma distrib
 rep,##<<Number of replicates
 obs.data##Formatted observed data
){
rain.data <- obs.data
  ls.monthly.simrain <- list()##List to store sim data for each month

  if(length(amount.param[1,]) == 1){ #If exponential

    for (i in 1:12){
      print(i)
      ls.monthly.simrain[[i]] <- data.frame(matrix(NA,nrow=length(rain.data[rain.data$month==i,2]),ncol=rep))
      for (j in 1:rep){
        # Create the occurrence binary series
        bin <- monthly_occurr_model(length(rain.data[rain.data$month==i,2]), occur.param[i,1], occur.param[i,2])
        # Declare monthly dataframe to store the SimRain
        pred <- rep(0,length(bin))
        # Loop to generate rainfall amount for rainny day
        for(k in which(bin == 1)){

          pred[k] <- rexp(1,amount.param[i,])
        }
        ls.monthly.simrain[[i]][,j] <- pred
      }
    }

    return(ls.monthly.simrain)

  }

  else if (length(amount.param[1,]) == 2){ #If gamma
    for (i in 1:12){
      print(i)
      ls.monthly.simrain[[i]] <- data.frame(matrix(NA,nrow=length(rain.data[rain.data$month==i,2]),ncol=rep))
      for (j in 1:rep){
        # Create the occurrence binary series
        bin <- monthly_occurr_model(length(rain.data[rain.data$month==i,2]), occur.param[i,1], occur.param[i,2])
        # Declare monthly dataframe to store the SimRain
        pred <- rep(0,length(bin))
        # Loop to generate rainfall amount for rainny day
        for(k in which(bin == 1)){

          pred[k] <- rgamma(1,shape=amount.param[i,1],rate=1/amount.param[i,2]) # rexp is the exponential distribution sampling function
        }
        ls.monthly.simrain[[i]][,j] <- pred
      }
    }
    return(ls.monthly.simrain)
  }
  else{
    print("check the input parameters data")
    return()
  }
}
