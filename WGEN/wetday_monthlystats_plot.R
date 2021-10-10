##Thien-10/06/2021
##Generate comparison plots for Mean, Std Dev, Skew monthly wet days
##Detail<<input are observed rainfall data frames and a list contains monthly simulated data frames
##Detail<<require additional functions: wet_day_stats() and pred.interval()
wetday_monthlystats_plot <- function
(monthly.obs, ##<< is a dataframe of observed rainfall at a gauged with the following format
 ##(column 1: Date; column 2: rainfall amount (mm); column 3: day; column 4:month; column 5: year)
 monthly.pred.list ##<<A list that contain 12 dataframes corresponding for 12 month
 ##each dataframe contains the simulated rainfall values for each month with length equal to observed data.frames
 ){
  
  #require(ggplot2)
  require(moments)
  
  ##Passing observed and simulated data from input
  rain.data <- monthly.obs
  ls.monthly.simrain <- monthly.pred.list
  ##Declare lists and data frame to be used to store the calculates statistic for each month and each replicates
  rep = length(ls.monthly.simrain[[1]][1,]) #Getting number of replicates
  monthstats.list <- list()## List that contains 12 dataframes of monthly statistics for each replicates
  obs.monthwetday.stats <- data.frame(matrix(NA,5,12)) ##This dataframe contains 5 monthly wet day stats for observed data
  colnames(obs.monthwetday.stats)<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")##Names the columns corresponding for each month
  rep.obs <- list()##List contains dataframes that is (rep) replicate of each rows of obs.monthwetday.stats for plotting purpose
  pred.stats <- list() ##List contains dataframes of each month and each replicates stats
  pred.CI.stats <- list()##LIst contains the 90% confidence limit and mean values of simulation sample
  
  ##Calculate simulated monthly wet day statistics 5 statistics: Mean, Skew, Std Dev, Min, Max for each month and for each rep
  for (i in 1:12){## Loop for 12 months
    df <- data.frame(matrix(NA,nrow = 5, ncol = rep)) #Declare property of the statistics dataframe: 
    monthstats.list[[i]] <- df #5 rows for 5 stats, rep columns for number of rep
    
    month.dat <- ls.monthly.simrain[[i]] ##Passing data of each month from the input list
    for(j in 1:rep){ ## Loop for number of replicates
      x<-rain.data[rain.data$month==i,] ## Getting for dates format to work
      x[,2]<-month.dat[,j] ## Replace the obs values with sim values
      monthstats.list[[i]][,j] <- wet_day_stats(x) ## calculating the statistics using the fn wet_day_stats()
    }
  }
  
  ##calculate observed monthly wet day statistics.
  for (i in 1:12){##Loop for 12 months
    obs.monthwetday.stats[,i] <- wet_day_stats(rain.data[rain.data$month==i,]) ##Calculate the statistics using fn wet_day_stats()
     
    }

  ##Gathering and preparing data for plotting
  for (i in 1:5){##Loop for each stats
    df1 <- data.frame(matrix(NA,rep,12))##declare property of dataframes (rep) number of rows and 12 columns
    pred.stats[[i]] <- df1 ##rep number of rows and 12 columns
    rep.obs[[i]] <- df1 ##replicates of obs stats with (rep) numbers of identical rows for 12 months
    #pred.CI.stats[[i]] <- data.frame(matrix(NA,3,12))##declare data frames contains upper,lower 90% CI and mean for each month
    for (j in 1:12){##Loop for each month
      pred.stats[[i]][,j] <- t(monthstats.list[[j]][i,])##Passing each stats for each month and each reps to pred.stats dataframe
      #pred.CI.stats[[i]][,j] <- CI90(pred.stats[[i]][,j])##Calculate 90% CI for the pred.stats dataframe using fn CI90()
    }
    colnames(pred.stats[[i]])<-c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec") ##Names the columns corresponding for each month 
    rep.obs[[i]] <- mapply(rep,obs.monthwetday.stats[i,],rep)##Replicate (rep) number of identical rows for obs stats
  }
  
  ##Start plotting
  title_plot1 <- c("Monthly wet day amounts (mm): means","Monthly wet day amounts (mm): skew","Monthly wet day amounts (mm): std dev")
  ylab_plot_2 <- c("mean montly wet day amounts (mm)","skew montly wet day amounts (mm)","std dev montly wet day amounts (mm)")
  
  for (i in 1:3){
    ##Type 1
    ##If statements are used to select the x and y limit of the plot
    if (min(pred.stats[[i]])<min(t(obs.monthwetday.stats[i,]))){
      min.range = min(pred.stats[[i]])
    } else {min.range = min(t(obs.monthwetday.stats[i,]))}
    if (max(pred.stats[[i]])>max(t(obs.monthwetday.stats[i,]))){
      max.range = max(pred.stats[[i]])
    } else {max.range = max(t(obs.monthwetday.stats[i,]))}
    ##Box plot function with modifiable whisker range 
    boxplot.ext(pred.stats[[i]],
                ylim=c(min.range,max.range),
                whiskersProb = c(0.05,0.95))
    ##Point the observed data on the plot
    points(t(obs.monthwetday.stats[i,]),col="red",pch=3)
    title(ylab=ylab_plot_2[i])
    
    ##Type 2: Error bar, using arrows() function to mimic error bar
    ## calculate the range of the arrow
    x<- sapply(pred.stats[[i]], pred.interval)
    ## Scatter plot the each 50%percentile point of the sim data to each point of the observed
    plot(x[2,],t(obs.monthwetday.stats[i,]),pch=3,col="blue",cex=0.5,
         xlim=c(min(x[1,]),max(x[3,])),
         xlab="simulated",ylab="observed")
    ##Plot the arrows with range of 90% prediction interval
    arrows(x0=x[1,],y0=t(obs.monthwetday.stats[i,]),x1=x[3,],y1=t(obs.monthwetday.stats[i,]),code=3, angle=90, length=0.0,col = "blue")     
    abline(coef = c(0,1),col="red")
    title(main=title_plot1[i])
  }

  ##ATTEMPT ON GGPLOT
#   month <- mapply(rep,c(1:12),rep)##Generate (rep) number of monthly strings c(1:2) for plotting purpose
#   ## Using ggplot2 to plot for 3 stats Mean, Skew and Std Dev
#   for (i in 1:3){##Loop for each stats
#     upper = mapply(rep,t(pred.CI.stats[[i]][3,]),rep)##Generate (rep) identical 90% upper CI for each stats
#     lower = mapply(rep,t(pred.CI.stats[[i]][1,]),rep)##Generate (rep) idenrical 90% lower CI for each stats
#     ## Using ggplot
#     plot_1 <- ggplot(data = NULL,aes(x = t(pred.stats[[i]]),y = t(rep.obs[[i]]))) + theme_bw()##start formatting plot_1
#     ls.names <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")## Data Point names for Plot_1
#     ##Set data names alternatively for odd and even months to obtain clearer presentation
#     for (k in seq(2,12,2)){ ##Loop for even months
#       plot_1 <- plot_1 + 
#         annotate(geom = "curve", x = t(pred.CI.stats[[i]][2,k]), y = t(rep.obs[[i]])[k,1], xend = t(pred.CI.stats[[i]][2,k]), yend = t(rep.obs[[i]])[k,1]-0.5, 
#                  curvature = 0, arrow = arrow(length = unit(0, "mm"))) +
#         annotate("text",x=t(pred.CI.stats[[i]][2,k]),y=t(rep.obs[[i]])[k,1]-0.5,label=ls.names[k],size=4)
#     }
#     for (k in seq(1,11,2)){##Loop for odd months
#       plot_1 <- plot_1 + 
#         annotate(geom = "curve", x = t(pred.CI.stats[[i]][2,k]), y = t(rep.obs[[i]])[k,1], xend = t(pred.CI.stats[[i]][2,k]), yend = t(rep.obs[[i]])[k,1]+0.5, 
#                  curvature = 0, arrow = arrow(length = unit(0, "mm"))) +
#         annotate("text",x=t(pred.CI.stats[[i]][2,k]),y=t(rep.obs[[i]])[k,1]+0.5,label=ls.names[k],size=4)
#     }
#     ##Continue format plot_1
#     plot_1 <- plot_1 + geom_abline(slope=1,intercept=0,color="red") +
#       labs(title = title_plot1[i],x = "simulated", y = "observed") +
#       theme(plot.title = element_text(hjust = 0.5)) + 
#       scale_x_continuous(n.breaks = 8) + scale_y_continuous(n.breaks = 8) +
#       expand_limits(x=c(min(t(rep.obs[[i]])),max(t(rep.obs[[i]]))),y=c(min(t(rep.obs[[i]])),max(t(rep.obs[[i]])))) +
#       geom_errorbar(aes(xmin=t(lower),xmax=t(upper)),width = 0.08,size = 1, color = "blue")
#     
#     print(plot_1)##Print plot_1
#     
#     ##Start formatting plot_2
#     plot_2 <-  ggplot(data=NULL,aes(x=month,y=t(pred.stats[[i]]))) + theme_bw() +
#       geom_errorbar(aes(ymin=upper,ymax=lower),width =0.2,size =0.8) +
#       geom_point(aes(x=c(1:12),y=t(obs.monthwetday.stats[i,]), color ="observed"),shape = 3,size = 1.5) +
#       geom_point(aes(x=c(1:12),y=t(pred.CI.stats[[i]][2,]), color = "simulated"),shape = 19,size = 1.9) +
#       labs(x = "month", y = ylab_plot_2[i]) +
#       scale_x_continuous(n.breaks = 12) + scale_y_continuous(n.breaks = 5) +
#       guides(colour=guide_legend(override.aes=list(shape=c(3,19)))) +
#       theme(legend.title = element_blank(), legend.position = c(0.9,0.9))+
#       scale_color_manual(
#         values=c("red", "blue"))
#     print(plot_2)##Print plot_2
#   }
#   
}