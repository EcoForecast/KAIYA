library(lubridate)
library(rjags)

siteid = 'HARV'
nee = targets_nee |> filter(site_id==siteid) |> filter(datetime<lubridate::date("2017-6-1"))
nee$datetime = lubridate::as_datetime(nee$datetime)

time = nee$datetime

plot(time, nee$observation, 'l')

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
                            n.iter = 10000)


time.rng = c(1,length(time))       ## adjust to zoom in and out
out <- as.matrix(jags.out)         ## convert from coda to matrix
x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(out[,x.cols],2,quantile,c(0.025,0.5,0.975))

plot(time,ci[2,],type='n',ylim=range(y,na.rm=TRUE),ylab="NEE",xlim=time[time.rng])
## adjust x-axis label to be monthly if zoomed
if(diff(time.rng) < 100){ 
  axis.Date(1, at=seq(time[time.rng[1]],time[time.rng[2]],by='month'), format = "%Y-%m")
}
ecoforecastR::ciEnvelope(time,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(time,y,pch="+",cex=0.2)

# Historical time-series fit -- linear model, HARV site

# we'll be using the ecoforecast tool, so let's install that, and load the rjags library in case we need it
library(rjags)
devtools::install_github("EcoForecast/ecoforecastR",force=TRUE)

# let's create a dataset that includes all of our variables that we can call it when fitting our model with JAGS
# this first run we'll be using the EFI temperature data, which has hourly resolution, while the NEE data has half hourly resolution. For now, we'll match each half hour measurement to the corresponding hourly temperature
# AT_half_hour <- list(parameter=rep(AT$parameter, each=2),datetime=as.POSIXlt(seq(from=as.POSIXct(AT$datetime[1]),to=as.POSIXct(AT$datetime[length(AT$datetime)]), by = "30 min")),variable=rep(AT$variable, each=2),prediction=rep(AT$prediction, each=2),site_id=rep(AT$site_id, each=2))
# let's match up the dates from the NEON temp data and the nee observations
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
