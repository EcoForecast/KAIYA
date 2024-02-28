library(neon4cast)
library(lubridate)
install.packages("rMR")
library(rMR)

## Data downloads

TODAYS_DATE = Sys.Date()

## LAI and FPAR

if(file.exists("01_LAI_FPAR_dwn.R"))      
  source("01_LAI_FPAR_dwn.R")

lai_download(TODAYS_DATE-30)
fpar_download(TODAYS_DATE-30)

## Temperature

if(file.exists("neontempdownload.R"))
  source("neontempdownload.R")

temp_data = temp_download(TODAYS_DATE)

# Save data to a file
temp_file = paste("NEON.Temp.","2018-01.",TODAYS_DATE,".HARV.RData", sep='')
file_path = file.path(getwd(), temp_file)

if(file.exists(file_path)){
}else{
save(subset,file=file_path)
}
