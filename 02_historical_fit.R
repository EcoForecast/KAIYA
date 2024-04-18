library(lubridate)
library(rjags)
library(dplyr)

# Set siteid
siteid = 'HARV'

# Get target variable (NEE)
if(file.exists('./Data/nee.RData')) {
  load('./Data/nee.RData')
} else if (file.exists("01_EFI_dwn.R")) {
  source("01_EFI_dwn.R")
  targets = download_targets()
  targets_nee = targets |> filter(variable=='nee') |> filter(site_id==siteid)
}

#### =============== Random walk =============== 
nee = targets_nee |> filter(site_id==siteid)
minDate = '2023-1-1'
maxDate = NULL
if(!is.null(minDate)) nee |> filter(datetime>lubridate::date(minDate))
if(!is.null(maxDate)) nee |> filter(datetime<lubridate::date(maxDate))
write.csv(nee, paste0('./Data/nee.', siteid, '.csv'))

nee$datetime = lubridate::as_datetime(nee$datetime)
time = nee$datetime
y = nee$observation
n_max = length(y)

RandomWalk = "
  model{

    #### Data Model
    for(t in 1:n){
      y[t] ~ dnorm(x[t],tau_obs)
    }

    #### Process Model
    for(t in 2:n){
      x[t]~dnorm(x[t-1],tau_add)
    }

    #### Priors
    x[1] ~ dnorm(x_ic,tau_ic)
    tau_obs ~ dgamma(a_obs,r_obs)
    tau_add ~ dgamma(a_add,r_add)
  }
  "
data <- list(y=y,n=length(y),      ## data
             x_ic=y[1],tau_ic=100,         ## initial condition prior
             a_obs=1,r_obs=1,                   ## obs error prior
             a_add=1,r_add=1                    ## process error prior
)

nchain = 3
init <- list()
for(i in 1:nchain){
  y.samp =sample(y,length(y),replace=TRUE)
  init[[i]] <- list(tau_add=1/var(diff(y.samp)),  ## initial guess on process precision
                    tau_obs=5/var(y.samp))        ## initial guess on obs precision
}

j.model   <- jags.model (file = textConnection(RandomWalk),
                         data = data,
                         inits = init,
                         n.chains = 3)

jags.out   <- coda.samples (model = j.model,
                            variable.names = c("tau_add","tau_obs"),
                            n.iter = 1000)
plot(jags.out)

jags.out   <- coda.samples (model = j.model,
                            variable.names = c("x","tau_add","tau_obs"),
                            n.iter = 5000)

time.rng = c(1,length(time))       ## adjust to zoom in and out
out <- as.matrix(jags.out)         ## convert from coda to matrix
x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(out[,x.cols],2,quantile,c(0.025,0.5,0.975))

# plot(time,ci[2,],type='n',ylim=range(y,na.rm=TRUE),ylab="NEE",xlim=time[time.rng])
# ## adjust x-axis label to be monthly if zoomed
# if(diff(time.rng) < 100){ 
#   axis.Date(1, at=seq(time[time.rng[1]],time[time.rng[2]],by='month'), format = "%Y-%m")
# }
# ecoforecastR::ciEnvelope(time,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
# points(time,y,pch="+",cex=0.2)
  
#### =============== Simple Dynamic Linear Model ===============
# if(file.exists(paste0('./Data/', 'HARV', '.met.historical.RData'))) {
#   load(paste0('./Data/', 'HARV', '.met.historical.RData'))
# } else if (file.exists("01_NOAA_dwn.R")) {
#   source("01_NOAA_dwn.R")
#   hist_met = weather_historical_download("2017-01-01", siteid)
# }
# 
# ## aggregation
# hist_met_aggregated <- hist_met %>% 
#   # mutate(datetime = lubridate::as_datetime(datetime)) %>% 
#   group_by(datetime, site_id, variable) |> 
#   summarize(prediction = mean(prediction),.groups = "drop") |> 
#   select(datetime, site_id, variable, prediction)
# 
# air_temp <- hist_met_aggregated |> filter(variable=='air_temperature')

####====================Amber's code====================

if(file.exists("01_NOAA_dwn.R"))
  source("01_NOAA_dwn.R")

temp_data = noaa_historical_download("HARV","air_temperature",Sys.Date()-365) # dated this way in order to get all data a year before today, this can be changed

## the following is leftover from when we were using EFI temp data
# temp_data = temp_data[c('mean_prediction', 'sd_prediction')] |> rename(
#   temp=mean_prediction,
#   temp_sd=sd_prediction
# )
# temp_data$datetime <- lubridate::as_datetime(temp_data$datetime)

# let's create a dataset that includes all of our variables that we can call it when fitting our model with JAGS
# this first run we'll be using the NOAA temperature data, which has hourly resolution, while the NEE data has half hourly resolution. For now, we'll match each half hour measurement to the corresponding hourly temperature
temp_half_hour <- list(datetime=as.POSIXlt(seq(from=as.POSIXct(temp_data$datetime[1]),to=as.POSIXct(temp_data$datetime[length(temp_data$datetime)]), by = "30 min")),temp=rep(temp_data$mean_prediction,each=2)[-length(temp_data$mean_prediction)],temp_sd=rep(temp_data$sd_prediction,each=2)[-length(temp_data$sd_prediction)],ensemble=rep(temp_data$ensemble, each=2)[-length(temp_data$ensemble)])

# Historical time-series fit -- linear model, HARV site

# let's match up the dates from the NEON temp data and the nee observations
temp_and_nee <- merge(temp_half_hour,nee,by='datetime')

# the following was a way to filter the data by year, if needed
# > filter(datetime>lubridate::date('2022-01-01'))
# just_2022 <- dplyr::filter(temp_and_nee,format(datetime, "%Y") == "2022")

data <- list(datetime=temp_and_nee$datetime,nee=temp_and_nee$observation,temp=temp_and_nee$temp,     ## data
             nee_ic=temp_and_nee$observation[1],tau_ic=100, ## initial condition prior
             a_obs=1,r_obs=1,           ## obs error prior
             a_add=1,r_add=1            ## process error prior
)
write.csv(temp_and_nee, './Data/combined.csv')

# now, let's use the ecoforecast package to run JAGS with a simple version of our model
nee.out <- ecoforecastR::fit_dlm(model=list(obs="nee",fixed="~ 1 + X + temp",n.iter=2000),data)
strsplit(nee.out$model,"\n",fixed = TRUE)[[1]] ## so we can verify what the JAGS code should look like
params <- window(nee.out$params,start=1000) ## remove burn-in
predicts <- window(nee.out$predict, start=1000)
plot(params)
summary(params)
cor(as.matrix(params))
pairs(as.matrix(params))

# now, we'll plot the data with the model and ci:
## confidence interval
out <- as.matrix(predicts)
write.csv(out[,(ncol(out)-32*48+1):ncol(out)], './Data/sim_lin_predicts.csv')
ci <- apply(out,2,quantile,c(0.025,0.5,0.975))
plot(temp_and_nee$datetime,ci[2,],type='n',xlab='date',ylim=c(-10,10),ylab="NEE")
ecoforecastR::ciEnvelope(temp_and_nee$datetime,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(temp_and_nee$datetime,temp_and_nee$observation,pch="+",cex=0.5)

out.params <- as.matrix(params)
write.csv(out.params, './Data/sim_lin_params.csv')

# a more complicated version:
# model <- ecoforecastR::ParseFixed("1 + AT + exp(-k*LAI*SW))
# where k,LAI,SW,AT forecasted?