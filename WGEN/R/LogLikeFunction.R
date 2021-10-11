gamma.loglike <- function(theta,obs.data){
  n = length(obs.data)
  loglike_i = vector(length = n)
  alpha = theta[1]; beta = theta[2]
  for (i in 1: n){
    loglike_i[i] = log(obs.data[i])
  }
  
  loglike_n = (alpha - 1) * sum(loglike_i) - 1/beta * sum(obs.data) - n * alpha *log(beta) - n * log(gamma(alpha))
  #-n * alpha * log(beta) - n * log(gamma(alpha)) - 1/beta * sum(obs.data) + (alpha - 1) * sum(loglike_i)
  if(is.nan(loglike_n)) {
    browser()
  }
  return(loglike_n)
}


expo.loglike <- function(theta,obs.data){
  n = length(obs.data)
  loglike_i = vector(length = n)
  
  for (i in 1:n){
    loglike_i[i] = obs.data[i]
  }
  loglike_n = n*log(theta) - theta*sum(loglike_i)
  if(is.nan(loglike_n)) {
    browser()
  }
  return(loglike_n)
}