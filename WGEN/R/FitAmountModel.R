#' Calibrate the amount model: a gamma or a exponential distribution
#'
#' @param obs.data observed daily TS of rainfall
#' @param mod User specify a gamma or a exponential distribution
#'
#' @return Amount model parameters
#'
#' @export

fitAmountModel<- function
(obs.data,##<<Formated observed data
 mod##<<Choosing which mod. Either: expo or gama
){

  if (mod == "expo"){

    amount.param <- data.frame(matrix(NA,nrow = 12, ncol = 1))

    for (i in 1:12){

      x1 <- obs.data[obs.data$month==i,] ##Passing monthly data to a temporary object

      x2<-subset(x1, x1$Value!=0) ##Subsetting day with 0 values

      ##<<Using optim() to do the MLE of model parameters
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

    return(amount.param)

  }
  else if(mod == "gama"){

    amount.param <- data.frame(matrix(NA,nrow = 12, ncol = 2))

    for (i in 1:12){

      x1 <- obs.data[obs.data$month==i,]

      x2<-subset(x1, x1$Value!=0)

      int.theta=c(0.1,0.1)    # Arbitrary initial values
      lower.theta=c(1e-5,1e-5)    # Lower bounds of parameters
      upper.theta=c(1,1e5)  # Upper bounds of parameters

      opt=optim(fn=gamma.loglike, # Function to be optimised
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

    return(amount.param)
  }
  else{
    print("Incorrect mod selection. Try again with either (expo) or (gama)")
    return()
  }
}
