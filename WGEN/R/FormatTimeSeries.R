#' Formatting a rainfall time series in to day, month, year
#'
#' @param x the Time Series
#'
#' @return A mutated Time Series with new columns for day, month and year
#'
#' @export
format_TimeSeries <- function(x){
  library(lubridate)
  library(dplyr)
  x <- x %>% #Set format for Date column
    mutate(Date_Time = as.Date(Date_Time, format = "%d/%m/%Y")) %>%
    mutate(month = month(Date_Time)) %>%
    mutate(year = year(Date_Time)) %>%
    mutate(day = day(Date_Time))
  return (x)
}
