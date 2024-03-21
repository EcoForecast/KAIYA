library(lubridate)
library(rjags)

# Set siteid
siteid = 'HARV'

# Get target variable (NEE)
if(file.exists('./Data/nee.RData')) {
  load('./Data/nee.RData')
} else if (file.exists("01_EFI_dwn.R")) {
  source("01_EFI_dwn.R")
  targets = download_targets()
  targets_nee = targets |> filter(variable=='nee' & site_id==siteid)
}

#### =============== Random walk ===============
nee = targets_nee |> filter(site_id==siteid) ### |> filter(datetime<lubridate::date("2017-6-1"))
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
                            variable.names = c("x","tau_add","tau_obs"),
                            n.iter = 1000)
plot(jags.out)

# jags.out   <- coda.samples (model = j.model,
#                             variable.names = c("x","tau_add","tau_obs"),
#                             n.iter = 10000)


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

# #### =============== Simple Dynamic Linear Model ===============
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
