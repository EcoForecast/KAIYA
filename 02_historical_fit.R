library(lubridate)
library(rjags)
library(dplyr)

# Set siteid
siteid = 'HARV'

## Get target variable (NEE)
if(file.exists('./Data/nee.RData')) {
  load('./Data/nee.RData')
} else if (file.exists("01_EFI_dwn.R")) {
  source("01_EFI_dwn.R")
  targets = download_targets()
  targets_nee = targets |> filter(variable=='nee') |> filter(site_id==siteid)
}

####====================Amber's code====================
if(file.exists('./Data/temp.hist.RData')) {
  load('./Data/temp.hist.RData')
} else if(file.exists("01_NOAA_dwn.R")){
  source("01_NOAA_dwn.R")
  temp_data = noaa_historical_download("HARV","air_temperature",Sys.Date()-365, Sys.Date())
}

semi_hourly_times = c()
mean_predictions = c()
sd_predictions = c()

for(i in 1:length(temp_data$datetime)){
  semi_hourly_times = c(semi_hourly_times, lubridate::as_datetime(temp_data$datetime[i]), lubridate::as_datetime(temp_data$datetime[i]+30*60))
  mean_predictions = c(mean_predictions, temp_data$mean_prediction[i], temp_data$mean_prediction[i])
  sd_predictions = c(sd_predictions, temp_data$sd_prediction[i], temp_data$sd_prediction[i])
}

temp_half_hour = data.frame(
  datetime = lubridate::as_datetime(semi_hourly_times),
  mean_prediction = mean_predictions,
  sd_prediction = sd_predictions
)

# let's create a dataset that includes all of our variables that we can call it when fitting our model with JAGS
# this first run we'll be using the NOAA temperature data, which has hourly resolution, while the NEE data has half hourly resolution. For now, we'll match each half hour measurement to the corresponding hourly temperature
# temp_half_hour <- list(datetime=as.POSIXlt(seq(from=as.POSIXct(temp_data$datetime[1]),to=as.POSIXct(temp_data$datetime[length(temp_data$datetime)]), by = "30 min")),temp=rep(temp_data$mean_prediction,each=2)[-length(temp_data$mean_prediction)],temp_sd=rep(temp_data$sd_prediction,each=2)[-length(temp_data$sd_prediction)],ensemble=rep(temp_data$ensemble, each=2)[-length(temp_data$ensemble)])

# Historical time-series fit -- linear model, HARV site

# let's match up the dates from the NEON temp data and the nee observations
# temp_and_nee <- merge(temp_half_hour,targets_nee,by='datetime')
combine_df <- merge(temp_half_hour,targets_nee,by='datetime')

# the following was a way to filter the data by year, if needed
# > filter(datetime>lubridate::date('2022-01-01'))
# just_2022 <- dplyr::filter(temp_and_nee,format(datetime, "%Y") == "2022")

data <- list(datetime=combine_df$datetime,nee=combine_df$observation,temp=combine_df$mean_prediction,     ## data
             nee_ic=combine_df$observation[1],tau_ic=100, ## initial condition prior
             a_obs=1,r_obs=1,           ## obs error prior
             a_add=1,r_add=1            ## process error prior
)
save(combine_df, file=file.path(getwd(), '/Data/combined.RData'))

# now, let's use the ecoforecast package to run JAGS with a simple version of our model
nee.out <- ecoforecastR::fit_dlm(model=list(obs="nee",fixed="~ 1 + X + temp",n.iter=10000),data)
strsplit(nee.out$model,"\n",fixed = TRUE)[[1]] ## so we can verify what the JAGS code should look like
params <- window(nee.out$params,start=5000) ## remove burn-in
predicts <- window(nee.out$predict, start=5000)
# plot(params)
summary(params)
cor(as.matrix(params))
# pairs(as.matrix(params))

# now, we'll plot the data with the model and ci:
## confidence interval
out <- as.matrix(predicts)
IC = out[,ncol(out)]
save(IC, file=file.path(getwd(), '/Data/sim_lin_IC.RData'))
ci <- apply(out,2,quantile,c(0.025,0.5,0.975))
filtered_ci = ci[,(ncol(ci)-2*48+1):ncol(ci)]
save(filtered_ci, file=file.path(getwd(), '/Data/sim_lin_ci.RData'))
plot(combine_df$datetime,ci[2,],type='n',xlab='date',ylim=c(-10,10),ylab="NEE")
ecoforecastR::ciEnvelope(combine_df$datetime,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(combine_df$datetime,combine_df$observation,pch="+",cex=0.5)

params <- as.matrix(params)
# write.csv(params, './Data/sim_lin_params.csv')
save(params, file=file.path(getwd(), '/Data/sim_lin_params.RData'))

# a more complicated version:
# model <- ecoforecastR::ParseFixed("1 + AT + exp(-k*LAI*SW))
# where k,LAI,SW,AT forecasted?