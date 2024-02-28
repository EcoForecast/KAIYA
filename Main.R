library(neon4cast)
library(lubridate)
install.packages("rMR")
library(rMR)

if(file.exists("01_LAI_FPAR_dwn.R"))      
  source("01_LAI_FPAR_dwn.R")

lai_download(Sys.Date()-30)
fpar_download(Sys.Date()-30)