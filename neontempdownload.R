install.packages("neonUtilities")
library(neonUtilities)

temp_download <- function(d) { #date input as YYYY-MM-DD, same as Sys.Date() gives
  temp_file <- paste("NEON.Temp.","2018-01.",substr(d,0,7),".HARV.RData", sep='')
  file_path <- file.path(getwd(), temp_file)
  if(file.exists(file_path)){print(file_path,'already exists')}
  else{
  download <- loadByProduct(dpID = "DP1.00003.001", site = "HARV", startdate = "2018-01", enddate = substr(d,0,7), tabl = "TAAT_30min", check.size = "F")
  dat <- download$TAAT_30min
  }
  return(dat)
}