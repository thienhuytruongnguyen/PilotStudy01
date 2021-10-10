pred.interval <- function(x,bound=c(0.05,0.95)){
  interval <- vector(length = 3)
  interval[1] <- quantile(x,prob = bound[1])
  interval[2] <- quantile(x,prob = 0.5)
  interval[3] <- quantile(x,prob = bound[2])
  return(interval)
}