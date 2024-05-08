library(neon4cast)
library(lubridate)
# install.packages("rMR")
# library(rMR)
library(rjags)

## Meta Data
team_info <- list(team_name = "KAIYA",
                  team_list = list(
                    list(individualName = list(givenName = "Katherine", 
                                               surName = "Rein"),
                         organizationName = "Boston University",
                         electronicMailAddress = "krein21@bu.edu"),
                    list(individualName = list(givenName = "Yuhe", 
                                               surName = "Chang"),
                         organizationName = "Boston University",
                         electronicMailAddress = "yhchang@bu.edu"),
                    list(individualName = list(givenName = "Ibbu", 
                                               surName = "Quraishi"),
                         organizationName = "Boston University",
                         electronicMailAddress = "quraiibr@bu.edu"),
                    list(individualName = list(givenName = "Amber", 
                                               surName = "Crenna-Armstrong"),
                         organizationName = "Boston University",
                         electronicMailAddress = "acrennaa@bu.edu"),
                    list(individualName = list(givenName = "Alex", 
                                               surName = "Coast"),
                         organizationName = "Boston University",
                         electronicMailAddress = "amocast@bu.edu"))
)

## Data downloads
TODAYS_DATE = Sys.Date()
siteid = 'HARV'

# ## LAI and FPAR
# if(file.exists("01_LAI_FPAR_dwn.R")){
#   source("01_LAI_FPAR_dwn.R")
#   lai_download(TODAYS_DATE-30)
#   fpar_download(TODAYS_DATE-30)
# } 
  
## EFI target variable
if(file.exists("01_EFI_dwn.R")){
  source("01_EFI_dwn.R")
  targets = download_targets()
  targets_nee = targets |> filter(variable=='nee') |> filter(site_id==siteid)
  targets_nee = targets_nee[order(targets_nee$datetime),]
  save(targets_nee, file = './Data/nee.RData')
  site_meta = download_site_meta()
  save(site_meta, file = './Data/site.RData')
  
  plot(targets_nee$datetime, targets_nee$observation, type='l', 
       main='NEE Time Series', xlab = 'Date', ylab = 'NEE')
}
  
## NOAA meteorological data
if(file.exists("01_NOAA_dwn.R")){
  source("01_NOAA_dwn.R")
  
  temp_data = noaa_historical_download(siteid,"air_temperature",TODAYS_DATE-365, TODAYS_DATE)
  temp_data = temp_data[order(temp_data$datetime),]
  save(temp_data, file='./Data/temp.hist.RData')
  
  plot(temp_data$datetime, temp_data$mean_prediction, 'l',
       main = 'Air Temperature', xlab = 'Date', ylab = 'Air temperature')
  
  temp_forecast = noaa_forecast_download(siteid,"air_temperature",
                                         as_date(max(targets_nee$datetime))+1)
  save(temp_data, file='./Data/temp.4cast.RData')
}
  
# download historical met data for HARV as a test

# siteid = "HARV"
# filename = paste0('./Data/', siteid, '.met.historical.RData')
# if(file.exists(filename)){} else {
#   hist_met = noaa_historical_download(siteid,Sys.Date()-80-365, Sys.Date()-80)
#   save(hist_met, file=filename)
# }
# 
# ## download historical met data for all sites
# # load('./Data/site.RData')
# # siteids = unique(site_meta$field_site_id)
# # for(siteid in siteids){
# #   ds = weather_historical_download("2024-03-01", siteid)
# #   filename = paste0('./Data/', siteid, '.met.historical.RData')
# #   save(ds, file=filename)
# # }
# 
# # update met data for HARV as a test
# siteid = "HARV"
# load(paste0('./Data/', siteid, '.met.historical.RData'))
# date_dwn = Sys.Date() # update data for today
# met = noaa_historical_download(siteid,,date_dwn-1,date_dwn)
# hist_met = rbind(hist_met, met)
# filename = paste0('./Data/', siteid, '.met.historical.RData')
# save(hist_met, file=filename)
# 
# ## time series plot of NEE for HARV
# load('./Data/nee.RData')
# nee = targets_nee |> filter(site_id=='HARV')
# plot(nee$datetime, nee$observation, type='l', 
#      main='NEE Time Series', xlab = 'Date', ylab = 'NEE')
# 
# ## time series plots of met covariates
# load(paste0('./Data/', 'HARV', '.met.historical.RData'))
# 
# # precipitation
# PRCP = hist_met[hist_met$variable == "precipitation_flux",]
# plot(PRCP$datetime, PRCP$prediction, 'l', 
#      main = 'Precipitaiton flux', xlab = 'Date', ylab = 'Precipitation')
# 
# # air temperature
# AT = hist_met[hist_met$variable == "air_temperature",]
# plot(AT$datetime, AT$prediction, 'l',
#      main = 'Air Temperature', xlab = 'Date', ylab = 'Air temperature')
# 
# # downwelling shortwave radiation
# SW = hist_met[hist_met$variable == "surface_downwelling_shortwave_flux_in_air",]
# plot(SW$datetime, SW$prediction, 'l', 
#      main = 'Downwelling Shortwave R', xlab = 'Date', ylab = 'Radiation')

# Historical fit
if(file.exists("02_historical_fit.R"))
  source("02_historical_fit.R")

plot(combine_df$datetime,ci[2,],type='n',xlab='date',ylim=c(-10,10),ylab="NEE")
ecoforecastR::ciEnvelope(combine_df$datetime,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(combine_df$datetime,combine_df$observation,pch="+",cex=0.5)

# Ensemble Forecast

if(file.exists("03_ensemble_4cast.R"))      
  source("03_ensemble_4cast.R")
plot.run()
N.IPDE.ci = apply(N.IPDE,2,quantile,c(0.025,0.5,0.975))
ecoforecastR::ciEnvelope(time2,N.IPDE.ci[1,],N.IPDE.ci[3,],col=col.alpha(N.cols[4],trans))
ecoforecastR::ciEnvelope(time2,N.IPD.ci[1,],N.IPD.ci[3,],col=col.alpha(N.cols[3],trans))
ecoforecastR::ciEnvelope(time2,N.IP.ci[1,],N.IP.ci[3,],col=col.alpha(N.cols[2],trans))
ecoforecastR::ciEnvelope(time2,N.I.ci[1,],N.I.ci[3,],col=col.alpha(N.cols[1],trans))
lines(time2,N.I.ci[2,],lwd=0.5)

# Model Validation

if(file.exists("O4_model_validation.R"))
  source("03_ensemble_4cast.R")
knitr::kable(stats)

plot(E,O,pch=".",xlab="ensemble",ylab='observed',main='NEE (umol/m2/sec)')
abline(0,1,col=2,lwd=2)
abline(NEE.ens.fit,col=3,lwd=3,lty=2)
legend("bottomright",legend=c('obs','1:1','reg'),col=1:3,lwd=3)

## Taylor diagrams: Ensemble
taylor.diagram(ref=O,model=E,normalize=TRUE,ref.sd=TRUE)
for(i in 1:nrow(NEE.ens)){
  taylor.diagram(ref=O,model=NEE.ens[i,qaqc],col=2,pch=".",add=TRUE,normalize=TRUE)
}
rlaplace = function(n,mu,b){
  return(mu + ifelse(rbinom(n,1,0.5),1,-1)*rexp(n,b))
}
beta = ifelse(O > 0,0.62+0.63*O,1.42-0.19*O) #Heteroskedasticity, parameters from Richardson et al 2006
for(i in 1:200){
  x = rlaplace(length(O),O,beta)
  taylor.diagram(ref=O,model=x,col=3,add=TRUE,normalize=TRUE)
}
legend("topright",legend=c("ens","obsUncert"),col=2:3,pch=20,cex=0.7)

plot(pval)   ## quantile 'residuals'
hist(pval,probability=TRUE) ## quantile distribution (should be flat)

mean(crps)
plot(crps)
hist(crps)

plot(time.of.day,crps.diurnal)
plot(O,crps)


# submission = false for now
source('./submit_forecast.R')
submit_forecast(result_df,team_info,submit=TRUE)

