library(ecoforecastR)
library(dplyr)
library(tidyverse)

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

if(file.exists('./Data/temp.4cast.RData')) {
  load('./Data/temp.4cast.RData')
} else if(file.exists("01_NOAA_dwn.R")) {
  source("01_NOAA_dwn.R")
  temp_forecast = noaa_forecast_download("HARV","air_temperature",as_date(max(combine_df$datetime))+1)
}

temp_forecast$datetime = as_datetime(temp_forecast$datetime)
n_ensemble = max(temp_forecast$parameter)+1

temp_ensemble <- matrix(NA, n_ensemble, (nrow(temp_forecast)/n_ensemble*2))
for(i in 1:n_ensemble){
  temp_df = temp_forecast[temp_forecast$parameter==(i-1),]
  ordered_df = temp_df[order(temp_df$datetime), ]
  temp_ensemble[i,] <- rep(ordered_df$prediction-273.15, each=2) # hourly to half-hourly
}
n_timesteps = ncol(temp_ensemble)

NT = 48
time = 1:(NT*2 + n_timesteps)    ## total time
time1 = 1:(NT*2)       ## calibration period
time2 = (NT*2+1):(NT*2 + n_timesteps)   ## forecast period
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
  N <- matrix(NA,n,n_timesteps)  ## storage
  Nprev <- IC           ## initialize
  for(t in 1:(n_timesteps)){
    mu = intercept + Nprev*(1+betax) + temperature[,t]*betatemp
    N[,t] <- rnorm(n,mu,Q)                         ## predict next step
    Nprev <- N[,t]                                  ## update IC
  }
  return(N)
}

## calculate mean of all inputs
temp.mean <- matrix(apply(temp_ensemble, 2, mean),1,(nrow(temp_forecast)/n_ensemble*2))
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

result_df <- data.frame(matrix(ncol = 9, nrow = 0))
predictions = c()
for(parameter in 1:100){
  for(i in 1:n_timesteps){
    result_df[nrow(result_df)+1,]=c('neon4cast','PT30','kaiya', 
                                    as.character(lubridate::as_datetime(min(temp_forecast$datetime)+30*60*(i-1))),
                                    as.character(combine_df$datetime[nrow(combine_df)]), 
                                    'HARV', 
                                    'ensemble',
                                    parameter,
                                    'nee'
                                    )
  }
  i_row = sample(1:Nmc, 1)
  predictions=c(predictions, N.IPDE[i_row,])
}
result_df = cbind(result_df, predictions)
col_names = c('project_id','duration','model_id', 'datetime', 'reference_datetime', 'site_id', 'family', 'parameter', 'variable', 'prediction')
colnames(result_df) <- col_names
