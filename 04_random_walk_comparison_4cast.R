library(lubridate)
library(rjags)
library(dplyr)
library(ecoforecastR)
library(dplyr)
library(tidyverse)

# Set siteid
TODAYS_DATE = Sys.Date()-365
siteid = 'HARV'

## Get target variable (NEE)
if(file.exists('./Data/nee.RData')) {
  load('./Data/nee.RData')
} else if (file.exists("01_EFI_dwn.R")) {
  source("01_EFI_dwn.R")
  targets = download_targets()
  targets_nee = targets |> filter(variable=='nee') |> filter(site_id==siteid)
}

combine_df <- targets_nee %>% filter(datetime >= TODAYS_DATE-365, datetime <= TODAYS_DATE)

data <- list(datetime=combine_df$datetime,nee=combine_df$observation,     ## data
             nee_ic=combine_df$observation[1],tau_ic=100, ## initial condition prior
             a_obs=1,r_obs=1,           ## obs error prior
             a_add=1,r_add=1            ## process error prior
)

# now, let's use the ecoforecast package to run JAGS with a simple version of our model
nee.out <- ecoforecastR::fit_dlm(model=list(obs="nee",fixed="~ 1 + X",n.iter=10000),data)
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
ci <- apply(out,2,quantile,c(0.025,0.5,0.975))
filtered_ci = ci[,(ncol(ci)-2*48+1):ncol(ci)]
plot(combine_df$datetime,ci[2,],type='n',xlab='date',ylim=c(-10,10),ylab="NEE")
ecoforecastR::ciEnvelope(combine_df$datetime,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(combine_df$datetime,combine_df$observation,pch="+",cex=0.5)
params <- as.matrix(params)
param.mean <- apply(params,2,mean)

#####################################################

## download temp data -- this is to make sure the dates match from our model's data validation
if(file.exists("01_NOAA_dwn.R")) {
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
semi_hourly_times = c()

for(i in 1:length(ordered_df$datetime)){
  semi_hourly_times = c(semi_hourly_times, lubridate::as_datetime(ordered_df$datetime[i]), lubridate::as_datetime(ordered_df$datetime[i]+30*60))
}

forecast_half_hour = data.frame(
  datetime = lubridate::as_datetime(semi_hourly_times)
)

##################################

date_forecast <- merge(targets_nee,forecast_half_hour,by='datetime',all.y=TRUE)$datetime
n_timesteps <- length(date_forecast)
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

forecastN <- function(IC, intercept, betax, Q=0, n=Nmc){
  N <- matrix(NA,n,n_timesteps)  ## storage
  Nprev <- IC           ## initialize
  for(t in 1:(n_timesteps)){
    mu = intercept + Nprev*(1+betax)
    N[,t] <- rnorm(n,mu,Q)                         ## predict next step
    Nprev <- N[,t]                                  ## update IC
  }
  return(N)
}

#### Deterministic Forecast
N.det <- forecastN(IC=mean(IC),
                   intercept=param.mean["betaIntercept"],
                   betax=param.mean["betaX"],
                   Q=0,  ## process error off
                   n=1)
lines(time2,N.det,col="purple",lwd=3)

#### Initial Condition uncertainty
prow = sample.int(nrow(params),Nmc,replace=TRUE)
N.I.rw <- forecastN(IC=IC[prow],
                 intercept=param.mean["betaIntercept"],
                 betax=param.mean["betaX"],
                 Q=0,  ## process error off
                 n=Nmc)
plot.run()
N.I.ci.rw = apply(N.I.rw,2,quantile,c(0.025,0.5,0.975))
ecoforecastR::ciEnvelope(time2,N.I.ci.rw[1,],N.I.ci.rw[3,],col=col.alpha(N.cols[1],trans))
lines(time2,N.I.ci.rw[2,],lwd=0.5)

save(N.I.rw,file=file.path(getwd(), '/Data/NI_randomwalk.RData'))

