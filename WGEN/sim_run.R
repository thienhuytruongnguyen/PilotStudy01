sim_run <- function
(occur.param,##<<Markov model parameters, PWW and PDW
 amount.param,##<<Amount model parameters, \lambda for exponential distrib or \alpha and \beta for gamma distrib
 rep=50##<<Number of replicates, default = 50
 ){
  
  ls.monthly.simrain <- list()##List to store sim data for each month
  
  if(length(amount.param[1,]) == 1){
    
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
  
  else if (length(amount.param[1,]) == 2){
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
