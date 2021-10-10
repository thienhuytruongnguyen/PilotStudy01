library(TimeSeriesManipulator)

#Import the observed data for two site
rain.data <- read.csv("41046_SILO_Rain.csv")
df1 <- readLines("month.txt") #List of monthly dataframe

#Formatting the obs data
rain.data <- format_TimeSeries(rain.data)

monthly.occur.param <- data.frame(matrix(nrow = 12,ncol = 4)) # This dataframe is to store the occurence parameters for each month
amount.param <- data.frame(matrix(NA,nrow = 12, ncol = 1))# This dataframe is to store the rainfall amount parameter Lambda
colnames(monthly.occur.param) <- c("PDW","PWW","PieW","PieD")

#CALIBRATE the occurrence and amount model parameters (PDW,PWW,Lambda)
for (i in 1:12){
  
  print(i)
  x <- occurr_stats(rain.data[rain.data$month==i,2]) #x is a temporary dataframe to store the outcome of function occurr_stats
  monthly.occur.param[i,] <- x[1,] #Store the occurrence parameters set for each month
  
  x1 <- rain.data[rain.data$month==i,]
  
  x2<-subset(x1, x1$Value!=0)
  
  int.theta=1    # Arbitrary initial values
  lower.theta=1e-5    # Lower bounds of parameters
  upper.theta=1e5  # Upper bounds of parameters

  opt=optim(fn=expo.loglike, # Function to be optimised
            par=int.theta,      # Initial Parameter Values
            obs=x2[,2],     # Synthetic data passed to obs argument in loglike.multinormal
            control=list(fnscale=-1,trace=2),  # Change from minimization to maximisation problem (fnscale), put some output
            lower=lower.theta,   # Lower bounds of parameters
            upper=upper.theta,    # Upper bounds of parameters
            method="L-BFGS-B"  # Quasi-Newton method that allows Box constraints see ?optim for details
  )
  est.theta=opt[["par"]] # Optimised parameter values
  amount.param[i,] <- est.theta # Store the optimise parameter in the predefined dataframe
  
}

# Simulate the monthly rainfall time series
for (i in 1:12){
  
  # Create the occurrence binary series
  bin <- monthly_occurr_model(length(rain.data[rain.data$month==i,2]), monthly.occur.param[i,1], monthly.occur.param[i,2])
  
  # Declare monthly dataframe to store the SimRain
  eval(parse(text = paste(df1[i],"<- data.frame(matrix(NA,nrow = length(bin), ncol = 1))")))
  
  # Loop to generate rainfall amount for rainny day
  for(j in 1: length(bin[,1])){
    
    # If statement to sample from the exponential distribution with calibrated Lambda in the previous section on a rainny day
    if(bin[j,] == 1){
      eval(parse(text = paste(
        df1[i],"[j,]","<- rexp(1,amount.param[i,])" # rexp is the exponential distribution sampling function
      )))
    }else{
      eval(parse(text = paste(
        df1[i],"[j,]","<- 0"
      )))
    }
  }
}

