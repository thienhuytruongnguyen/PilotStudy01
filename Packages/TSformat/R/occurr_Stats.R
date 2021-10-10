#' Calculate the monthly occurrence statistic of a rainfall Time Series
#'
#' @param monthlydata The monthly Time Serie
#'
#' @return A dataframe containing 4 values: Probability of Dry-Wet Event "P_DW", Probability of Wet-Wet Event "PWW", Marginal Probability of Wet Days "Pie_W, Marginal Probability of Dry Days "Pie_D"
#'
#' @export
occurr_stats <- function(monthlydata){

  #Storing only rain data to a data frame for processing
  rain.only.month <- data.frame(matrix(NA,nrow = length(monthlydata[,2])))
  rain.only.month <- monthlydata[,2]

  P_dw = 0
  P_ww = 0
  Pie_w = 0
  Pie_d = 0
  P_c <- data.frame(matrix(NA,nrow = (length(rain.only.month)-1)))
  #Loop for each site

    #print(s)
    n_dw = 0 #Number of dw events at site
    n_ww = 0 #Number of ww events at site
    n_dd = 0 #Number of dd events at site
    n_wd = 0 #Number of wd events at site
    for (i in 2:length(rain.only.month)){
      if (rain.only.month[i] > 0.1 && rain.only.month[i-1] <= 0.1){
        n_dw = n_dw + 1
      }
      if (rain.only.month[i] > 0.1 && rain.only.month[i-1] > 0.1){
        n_ww = n_ww + 1
      }
      if (rain.only.month[i] <= 0.1 && rain.only.month[i-1] <= 0.1){
        n_dd = n_dd + 1
      }
      if (rain.only.month[i] <= 0.1 && rain.only.month[i-1] > 0.1){
        n_wd = n_wd + 1
      }

      #Calculate the conditional probability
      P_dw = n_dw/(n_dw+n_dd)
      P_ww = n_ww/(n_ww+n_wd)
      #Unconditional probabilty
      Pie_w = P_dw/(1 + P_dw - P_ww)
      Pie_d = 1 - Pie_w
    }
    sumary <- data.frame(matrix(NA,nrow = 1, ncol = 4))
    sumary[,1] = P_dw ; sumary[,2] = P_ww ; sumary[,3] = Pie_w ; sumary[,4] = Pie_d
    colnames(sumary) <- c("P_DW","P_WW","Pie_W","Pie_D")
       return(sumary)
}
