library(TimeSeriesManipulator)

#------------------------------Import the observed data for two site---------------------------------#
rain.data <- read.csv("Rainfall.csv")
#----------------------------------------------------------------------------------------------------#


#-------------------------------------Formatting the obs data----------------------------------------#
rain.data <- format_TimeSeries(rain.data)
#----------------------------------------------------------------------------------------------------#


#---------------------------------------Declaring objects--------------------------------------------#
occur.param <- data.frame() # This dataframe is to store the occurrence parameters for each month
amount.param <- data.frame()# This dataframe is to store the rainfall amount parameter Lambda
#----------------------------------------------------------------------------------------------------#


#--------------CALIBRATE the occurrence and amount model parameters (PDW,PWW,Lambda)-----------------#
occur.param <- calibrate.markov_model(obs.data = rain.data)
amount.param <- calibrate.amount_model(obs.data = rain.data,mod = "gama")
#----------------------------------------------------------------------------------------------------#


#---------------Simulate the monthly rainfall time series--------------------------------------------#
ls.month.sim <- sim_run(occur.param = occur.param, amount.param = amount.param, rep = 100)
#----------------------------------------------------------------------------------------------------#


#---------------------------------Write the parameters to .csv file----------------------------------#
write.csv(amount.param,file = "amount_para.csv")
write.csv(amount.param.gamma, file = "amount_para_gama.csv")
write.csv(monthly.occur.param, file = "tranprob.csv")
#----------------------------------------------------------------------------------------------------#


#-----------------------------Plotting---------------------------------------------------------------#
wetday_monthlystats_plot(rain.data, ls.month.sim) #Plot monthly wet day statistics (mean, skew, std dev)
monthlytotal_stats_plot(rain.data, ls.month.sim)
