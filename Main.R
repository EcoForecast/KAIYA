library(neon4cast)
library(lubridate)
install.packages("rMR")
library(rMR)
library(dplyr)

## Data downloads

TODAYS_DATE = Sys.Date()

## LAI and FPAR
## time resolution of 

if(file.exists("01_LAI_FPAR_dwn.R")){source("01_LAI_FPAR_dwn.R")}      

lai = lai_download(TODAYS_DATE-30)
fpar = fpar_download(TODAYS_DATE-30)

## Temperature
## time resolution of 1/2 hour

if(file.exists("neontempdownload.R")){
  source("neontempdownload.R")}

temp_data = temp_download(TODAYS_DATE)

# Save data to a file
temp_file = paste("NEON.Temp.","2018-01.",substr(TODAYS_DATE,0,7),".HARV.RData", sep='')
temp_path = file.path(getwd(), 'Data',temp_file)

if(file.exists(temp_path)){}else{
save(temp_data,file=temp_path)
}

## EFI target variable
## time resolution of 1/2 hour

if(file.exists("01_EFI_dwn.R")){
  source("01_EFI_dwn.R")}

targets = download_targets()
targets_nee = targets |> filter(variable=='nee')
save(targets_nee, file = './Data/nee.RData')

site_meta = download_site_meta()
save(site_meta, file = './Data/site.RData')

## NOAA meteorological data
## time resolution of hour

if(file.exists("01_NOAA_dwn.R")){
  source("01_NOAA_dwn.R")}

# download historical met data for HARV as a test

siteid = "HARV"
filename = paste0('./Data/', siteid, '.met.historical.RData')
if(file.exists(filename)){} else {
  hist_met = weather_historical_download("2021-01-01", siteid)
  save(hist_met, file=filename)
}

## download historical met data for all sites
# load('./Data/site.RData')
# siteids = unique(site_meta$field_site_id)
# for(siteid in siteids){
#   ds = weather_historical_download("2024-03-01", siteid)
#   filename = paste0('./Data/', siteid, '.met.historical.RData')
#   save(ds, file=filename)
# }

# update met data for HARV as a test
siteid = "HARV"
load(paste0('./Data/', siteid, '.met.historical.RData'))
date_dwn = Sys.Date() # update data for today
met = weather_download(date_dwn, siteid)
hist_met = rbind(hist_met, met)
filename = paste0('./Data/', siteid, '.met.historical.RData')
save(hist_met, file=filename)

## load data, define variables, time series plots

## time series plot of NEE for HARV
load('./Data/nee.RData')
nee = targets_nee |> filter(site_id=='HARV')
plot(nee$datetime, nee$observation, type='l', 
     main='NEE Time Series', xlab = 'Date', ylab = 'NEE')

## time series plots of met covariates
load(paste0('./Data/', 'HARV', '.met.historical.RData'))

# precipitation
PRCP = hist_met[hist_met$variable == "precipitation_flux",]
plot(PRCP$datetime, PRCP$prediction, 'l', 
     main = 'Precipitation flux', xlab = 'Date', ylab = 'Precipitation')

# air temperature
AT = hist_met[hist_met$variable == "air_temperature",]
plot(AT$datetime, AT$prediction, 'l',
     main = 'Air Temperature', xlab = 'Date', ylab = 'Air temperature')

# downwelling shortwave radiation
SW = hist_met[hist_met$variable == "surface_downwelling_shortwave_flux_in_air",]
plot(SW$datetime, SW$prediction, 'l', 
     main = 'Downwelling Shortwave R', xlab = 'Date', ylab = 'Radiation')

# NEON temperature
load('./Data/NEON.Temp.2018-01.2024-03.HARV.RData')
temp = list(datetime = temp_data$startDateTime,temp = temp_data$tempTripleMean)

# Historical time-series fit -- HARV site

# we'll be using the ecoforecast tool, so let's install that, and load the rjags library in case we need it
library(rjags)
devtools::install_github("EcoForecast/ecoforecastR",force=TRUE)

# let's create a dataset that includes all of our variables that we can call it when fitting our model with JAGS
# this first run we'll be using the EFI temperature data, which has hourly resolution, while the NEE data has half hourly resolution. For now, we'll match each half hour measurement to the corresponding hourly temperature
# AT_half_hour <- list(parameter=rep(AT$parameter, each=2),datetime=as.POSIXlt(seq(from=as.POSIXct(AT$datetime[1]),to=as.POSIXct(AT$datetime[length(AT$datetime)]), by = "30 min")),variable=rep(AT$variable, each=2),prediction=rep(AT$prediction, each=2),site_id=rep(AT$site_id, each=2))
#  let's match up the dates from the NEON temp data and the nee observations
temp_and_nee <- merge(temp,nee,by='datetime')
just_2022 <- dplyr::filter(temp_and_nee,format(datetime, "%Y") == "2022")
data <- list(datetime=just_2022$datetime,nee=just_2022$observation,temp=just_2022$temp,     ## data
             nee_ic=just_2022$observation[1],tau_ic=100, ## initial condition prior
             a_obs=1,r_obs=1,           ## obs error prior
             a_add=1,r_add=1            ## process error prior
)

# now, let's use the ecoforecast package to run JAGS with a simple version of our model
nee.out <- ecoforecastR::fit_dlm(model=list(obs="nee",fixed="~ 1 + X + temp"),data)
params <- window(nee.out$params,start=1000) ## remove burn-in
plot(params)
summary(params)
cor(as.matrix(params))
pairs(as.matrix(params))

# now, we'll plot the data with the model and ci:
## confidence interval
out <- as.matrix(nee.out$predict)
ci <- apply(out,2,quantile,c(0.025,0.5,0.975))
plot(just_2022$datetime,ci[2,],type='n',xlab='date',ylim=c(-10,10),ylab="NEE")
ecoforecastR::ciEnvelope(just_2022$datetime,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(just_2022$datetime,just_2022$observation,pch="+",cex=0.5)

# a more complicated version:
# model <- ecoforecastR::ParseFixed("1 + AT + exp(-k*LAI*SW))
# where we predict nee and k?
