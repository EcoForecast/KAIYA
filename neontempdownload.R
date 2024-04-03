install.packages("neonUtilities")
library(neonUtilities)

temp_download <- function(d) { #date input as YYYY-MM-DD, same as Sys.Date() gives
  download <- loadByProduct(dpID = "DP1.00003.001", site = "HARV", startdate = "2018-01", enddate = substr(Sys.Date(),0,7), tabl = "TAAT_30min", check.size = "F")
  dat <- download$TAAT_30min
  return(dat)
}
