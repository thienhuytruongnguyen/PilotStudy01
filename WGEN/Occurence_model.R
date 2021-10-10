#-----------------------------------------------------------------------------------------------------------------------#
#Created by Thien 03/06/2021
#This is a Markov Chain First order model to generate a binary series for monthly rainfall occurrence
#The output is in the form of a dataframe contain the binary series of monthly rainfall occurrence
#The inputs needed are the length of the monthly observed rainfall (N) and the transitional probabilities PDW and PWW
#-----------------------------------------------------------------------------------------------------------------------#

monthly_occurr_model <- function(N,PDW,PWW){
  ## Assume P_c = PDW for January
  P_c = PDW
  
  #Assign an dataframe to store the occurrence binary series
  x <- vector(length = N)
  
  #generate a uniform random series U[0,1] to force the occurrence binary series
  U_t <- runif(N,0,1)
  
  for (j in 1:length(U_t)){ #loop for generating the binary occurrence time series
    if (U_t[j] < P_c){ #If statement to sample a 1 or 0 rainfall occurrence from the uniform random series U_t and P_C
      x[j] = 1
    }else{x[j] = 0}
    
    if(x[j] == 1){ #If statement to Update the P_c for the next sampling
      P_c = PWW
    }else{P_c = PDW}
  }
  return(x)
}
