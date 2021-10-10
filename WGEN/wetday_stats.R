wet_day_stats <- function(x){
  
  require(moments)
  pred.stats <- vector(length = 5)
  
  wet.day.dat <- subset(x,x[,2]!=0)
  
  pred.stats[1] = mean(wet.day.dat[,2])
  pred.stats[2] = skewness(wet.day.dat[,2])
  pred.stats[3] = sd(wet.day.dat[,2])
  pred.stats[4] = min(wet.day.dat[,2])
  pred.stats[5] = max(wet.day.dat[,2])
  
  return(pred.stats)
}