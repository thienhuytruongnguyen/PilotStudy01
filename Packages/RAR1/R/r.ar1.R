#' The auto regressive lag 1 model
#'
#' @param n The sample length
#' @param mu The mean
#' @param sigma The Standard Deviation
#' @param phi The Correlation coefficient
#'
#' @return A sample of random variables
#'
#' @examples
#' r.ar1(n = 10, mu = 0, sigma = 1, phi = 0.1)
#'
#' @export
r.ar1=function(n,mu,sigma,phi){

  # Set-up vector
  pred=vector(length=n)

  # Use marginal distibution to sample first point
  marg.sd=sigma/sqrt(1-phi^2)
  pred[1]=rnorm(n=1,mean=mu,sd=marg.sd)

  #Sample eta's as vector
  eta=rnorm(n=n,mean=0,sd=sigma)

  # Produce ar(1) predictions
  for (i in 2:n){
    pred[i]=mu+phi*(pred[i-1]-mu)+eta[i]
  }

  return(pred)

}
