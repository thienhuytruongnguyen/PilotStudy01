CI90 <- function(x){
  CI <- vector(length=3)
  z = 1.64
  CI[1] = mean(x) - 1.64*(sd(x)/sqrt(length(x)))
  CI[3] = mean(x) + 1.64*(sd(x)/sqrt(length(x)))
  CI[2] = mean(x)
  return(CI)
}

