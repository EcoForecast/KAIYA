library(ecoforecastR)
library(dplyr)
library(tidyverse)

# if(file.exists('./Data/combined.csv')) {
#   combine_df = read.csv('./Data/combined.csv')[-1]
# }
# if(file.exists('./Data/sim_lin_params.csv')) {
#   params = read.csv('./Data/sim_lin_params.csv')[-1]
# }

if(file.exists('./Data/combined.RData')) {
  load('./Data/combined.RData')
}
if(file.exists('./Data/sim_lin_params.RData')) {
  load('./Data/sim_lin_params.RData')
}
if(file.exists('./Data/sim_lin_ci.RData')) {
  load('./Data/sim_lin_ci.RData')
  }

#### Download temperature forecast
if(file.exists("01_NOAA_dwn.R"))
  source("01_NOAA_dwn.R")

temp_forecast = noaa_forecast_download("HARV","air_temperature",Sys.Date()-80)
temp_forecast$datetime = as_datetime(temp_forecast$datetime)

#### Convert hourly to half-hourly
semi_hourly_times = c()
mean_predictions = c()
sd_predictions = c()
for(i in 1:length(temp_forecast$datetime)){
  semi_hourly_times = c(semi_hourly_times, lubridate::as_datetime(temp_forecast$datetime[i]), lubridate::as_datetime(temp_forecast$datetime[i]+30*60))
  mean_predictions = c(mean_predictions, temp_forecast$mean_prediction[i], temp_forecast$mean_prediction[i])
  sd_predictions = c(sd_predictions, temp_forecast$sd_prediction[i], temp_forecast$sd_prediction[i])
}

temp_forecast = data.frame(
  datetime = lubridate::as_datetime(semi_hourly_times),
  mean_prediction = mean_predictions,
  sd_prediction = sd_predictions
)

NT = 48
time = 1:(NT*2 + nrow(temp_forecast))    ## total time
time1 = 1:(NT*2)       ## calibration period
time2 = (NT*2+1):(NT*2 + nrow(temp_forecast))   ## forecast period
ylim=c(-20,20)
Nmc = 1000         ## set number of Monte Carlo draws
N.cols <- c("black","red","green","blue","orange") ## set colors
trans <- 0.5

plot.run <- function(){
  plot(time,time,type='n',ylim=ylim,ylab="N")
  ecoforecastR::ciEnvelope(time1,filtered_ci[1,],filtered_ci[3,],col=col.alpha("lightBlue",0.6))
  lines(time1,filtered_ci[2,],col="blue")
  # points(time1,last2days[s,])
}
plot.run()

forecastN <- function(IC, temperature, intercept, betax, betatemp, Q=0, n=Nmc){
  N <- matrix(NA,n,(nrow(temp_forecast)))  ## storage
  Nprev <- IC           ## initialize
  for(t in 1:(nrow(temp_forecast))){
    mu = intercept + Nprev*(1+betax) + temperature[,t]*betatemp
    N[,t] <- rnorm(n,mu,Q)                         ## predict next step
    Nprev <- N[,t]                                  ## update IC
  }
  return(N)
}

## calculate mean of all inputs
temp_ensemble <- matrix(NA, 2000, (nrow(temp_forecast)))
for(i in 1:(nrow(temp_forecast))){
  temp_ensemble[,i] <- rnorm(50, temp_forecast$mean_prediction[i], temp_forecast$sd_prediction[i])
}
temp.mean <- matrix(apply(temp_ensemble, 2, mean),1,(nrow(temp_forecast)))
## parameters
param.mean <- apply(params,2,mean)
load('./Data/sim_lin_IC.RData')

#### Deterministic Forecast
N.det <- forecastN(IC=mean(IC),
                   temperature=temp.mean,
                   intercept=param.mean["betaIntercept"],
                   betax=param.mean["betaX"],
                   betatemp=param.mean["betatemp"],
                   Q=0,  ## process error off
                   n=1)
lines(time2,N.det,col="purple",lwd=3)

#### Initial Condition uncertainty
prow = sample.int(nrow(params),Nmc,replace=TRUE)
N.I <- forecastN(IC=IC[prow],
                 temperature=temp.mean,
                 intercept=param.mean["betaIntercept"],
                 betax=param.mean["betaX"],
                 betatemp=param.mean["betatemp"],
                 Q=0,  ## process error off
                 n=Nmc)
plot.run()
N.I.ci = apply(N.I,2,quantile,c(0.025,0.5,0.975))
ecoforecastR::ciEnvelope(time2,N.I.ci[1,],N.I.ci[3,],col=col.alpha(N.cols[1],trans))
lines(time2,N.I.ci[2,],lwd=0.5)

#### Parameter uncertainty
N.IP <- forecastN(IC=IC[prow],
                  temperature=temp.mean,
                  intercept=params[prow, "betaIntercept"],
                  betax=params[prow,"betaX"],
                  betatemp=params[prow,"betatemp"],
                  Q=0,  ## process error off
                  n=Nmc)
## Plot run
plot.run()
N.IP.ci = apply(N.IP,2,quantile,c(0.025,0.5,0.975))
ecoforecastR::ciEnvelope(time2,N.IP.ci[1,],N.IP.ci[3,],col=col.alpha(N.cols[2],trans))
ecoforecastR::ciEnvelope(time2,N.I.ci[1,],N.I.ci[3,],col=col.alpha(N.cols[1],trans))
lines(time2,N.I.ci[2,],lwd=0.5)

#### Driver uncertainty
drow = sample.int(nrow(temp_ensemble),Nmc,replace=TRUE)
N.IPD <- forecastN(IC=IC[prow],
                   temperature=temp_ensemble[drow,],
                   intercept=params[prow, "betaIntercept"],
                   betax=params[prow,"betaX"],
                   betatemp=params[prow,"betatemp"],
                   Q=0,  ## process error off
                   n=Nmc)
## Plot run
plot.run()
N.IPD.ci = apply(N.IPD,2,quantile,c(0.025,0.5,0.975))
ecoforecastR::ciEnvelope(time2,N.IPD.ci[1,],N.IPD.ci[3,],col=col.alpha(N.cols[3],trans))
ecoforecastR::ciEnvelope(time2,N.IP.ci[1,],N.IP.ci[3,],col=col.alpha(N.cols[2],trans))
ecoforecastR::ciEnvelope(time2,N.I.ci[1,],N.I.ci[3,],col=col.alpha(N.cols[1],trans))
lines(time2,N.I.ci[2,],lwd=0.5)

#### process error samples
Qmc <- 1/sqrt(params[prow,"tau_add"])  ## convert from precision to standard deviation
N.IPDE <- forecastN(IC=IC[prow],
                    temperature=temp_ensemble[drow,],
                    intercept=params[prow, "betaIntercept"],
                    betax=params[prow,"betaX"],
                    betatemp=params[prow,"betatemp"],
                    Q=Qmc,
                    n=Nmc)
## Plot run
plot.run()
N.IPDE.ci = apply(N.IPDE,2,quantile,c(0.025,0.5,0.975))
ecoforecastR::ciEnvelope(time2,N.IPDE.ci[1,],N.IPDE.ci[3,],col=col.alpha(N.cols[4],trans))
ecoforecastR::ciEnvelope(time2,N.IPD.ci[1,],N.IPD.ci[3,],col=col.alpha(N.cols[3],trans))
ecoforecastR::ciEnvelope(time2,N.IP.ci[1,],N.IP.ci[3,],col=col.alpha(N.cols[2],trans))
ecoforecastR::ciEnvelope(time2,N.I.ci[1,],N.I.ci[3,],col=col.alpha(N.cols[1],trans))
lines(time2,N.I.ci[2,],lwd=0.5)

result_df <- data.frame(matrix(ncol = 10, nrow = 0))
col_names = c('project_id','duration','model_id', 'datetime', 'reference_datetime', 'site_id', 'family', 'parameter', 'variable', 'prediction')
colnames(result_df) <- col_names

for(i in 1:(nrow(temp_forecast))){
  for(parameter in 1:30){
    result_df[nrow(result_df)+1,]=c('neon4cast','PT30','kaiya', 
                                    as.character(temp_forecast$datetime[i]), 
                                    as.character(combine_df$datetime[nrow(combine_df)]), 
                                    'HARV', 
                                    'ensemble',
                                    parameter,
                                    'nee',
                                    sample(N.IPDE[,i], 1))
  }
}

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

# submission = false for now
source('./submit_forecast.R')
submit_forecast(result_df,team_info,submit=TRUE)
