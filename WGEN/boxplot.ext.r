boxplot.ext=function(x,whiskersProb=c(0.025,0.975),...){
# Draws a boxplot with the whiskers at the probability limits, provided by whiskersProb
# x can be vector or data.frame
# whiskers Prob values are probabilities for the lower and upper whisker (respectively)


x.stats=boxplot(x,plot=FALSE)
if ( is.data.frame(x) | is.matrix(x)) {
  y=x
  for (j in 1:ncol(x)) {
    y=sort(x[,j])
    if (length(na.omit(y))!=0) {
      x.stats$stats[1,j]=quantile(y,prob=whiskersProb[1])
      x.stats$stats[5,j]=quantile(y,prob=whiskersProb[2])
    }
  }
  }
else if (is.vector(x)) { # Assume x is a vector
  y=sort(x)
  if (length(na.omit(y))!=0) {
    x.stats$stats[1,1]=quantile(y,prob=whiskersProb[1])
    x.stats$stats[5,1]=quantile(y,prob=whiskersProb[2])
    }
  }
else
{
  print("type of x is not supported in boxplot.ext")
  return()
}

bxp(z=x.stats,...)

}


