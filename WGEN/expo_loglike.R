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