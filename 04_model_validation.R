## load libraries
install.packages("plotrix",repos = "http://cran.us.r-project.org")
library("plotrix")
library(rpart)
library(scoringRules)
library(ecoforecastR)
library(dplyr)
library(tidyverse)

## Data downloads
TODAYS_DATE = Sys.Date()-365
siteid = 'HARV'

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

## source historical fit

if(file.exists("02_historical_fit.R"))
  source("02_historical_fit.R")

## run forecast

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


### Model validation

# match dates from data and results
# NEE.ens.bar <- result_df %>% select(parameter,prediction,datetime) %>% mutate(ensemble = as.numeric(parameter),nee = as.numeric(prediction)) %>% group_by(datetime) %>% summarize(nee = -mean(nee))
NEE.ens <- N.I
NEE.ens.bar <- N.I.ci[2,]
# nee_data <- targets_nee %>% filter(datetime >= min(ordered_df$datetime) & datetime <= max(ordered_df$datetime))

## so, it turns out that there are some data missing from the targets_nee dataset
semi_hourly_times = c()

for(i in 1:length(ordered_df$datetime)){
  semi_hourly_times = c(semi_hourly_times, lubridate::as_datetime(ordered_df$datetime[i]), lubridate::as_datetime(ordered_df$datetime[i]+30*60))
}

forecast_half_hour = data.frame(
  datetime = lubridate::as_datetime(semi_hourly_times)
)
nee_data <- merge(targets_nee,forecast_half_hour,by='datetime',all.y=TRUE)

## Calculate ensemble means & apply QAQC
qaqc <- complete.cases(nee_data)
E <- NEE.ens.bar[qaqc]
O <- nee_data$observation[qaqc]

## Model vs obs regressions
NEE.ens.fit = lm(O ~ E)

## performance stats
stats = as.data.frame(matrix(NA,4))
rownames(stats) <- c("RMSE","Bias","cor","slope")
colnames(stats) <- c("ens")
stats["RMSE",'ens'] = sqrt(mean((E-O)^2))
stats['Bias','ens'] = mean(E-O)
stats['cor','ens']  = cor(E,O)
stats['slope','ens'] = coef(NEE.ens.fit)[2]
knitr::kable(stats)

## predicted-observed
plot(E,O,pch=".",xlab="ensemble",ylab='observed',main='NEE (umol/m2/sec)')
abline(0,1,col=2,lwd=2)
abline(NEE.ens.fit,col=3,lwd=3,lty=2)
legend("bottomright",legend=c('obs','1:1','reg'),col=1:3,lwd=3)

## Taylor diagrams: Ensemble

taylor.diagram(ref=O,model=E,normalize=TRUE,ref.sd=TRUE)
for(i in 1:nrow(NEE.ens)){
  taylor.diagram(ref=O,model=NEE.ens[i,qaqc],col=2,pch=".",add=TRUE,normalize=TRUE)
}

## add data uncertainty
rlaplace = function(n,mu,b){
  return(mu + ifelse(rbinom(n,1,0.5),1,-1)*rexp(n,b))
}
beta = ifelse(O > 0,0.62+0.63*O,1.42-0.19*O) #Heteroskedasticity, parameters from Richardson et al 2006
for(i in 1:200){
  x = rlaplace(length(O),O,beta)
  taylor.diagram(ref=O,model=x,col=3,add=TRUE,normalize=TRUE)
}
legend("topright",legend=c("ens","obsUncert"),col=2:3,pch=20,cex=0.7)


# Bayesian p-values
O <- nee_data$observation
pval = 0
for(i in 1:ncol(NEE.ens)){
  pval[i] = sum(O[i] > -NEE.ens[,i])/nrow(NEE.ens)  ## quantiles of the particle filter
}
plot(pval)   ## quantile 'residuals'
hist(pval,probability=TRUE) ## quantile distribution (should be flat)

### CRPS
# note: I attempted to run this code without performing qaqc, but since the data has missing values, I had to do the qaqc indexed matrices
O <- O[qaqc]
crps = scoringRules::crps_sample(y = O, dat = t(NEE.ens[,qaqc]))
mean(crps)
plot(crps)
hist(crps)


### diurnal cycle of model skill
crps.diurnal = tapply(crps,lubridate::minute(nee_data$datetime[qaqc]),mean)
time.of.day = as.numeric(names(crps.diurnal))
plot(time.of.day,crps.diurnal)

## CRPS skill as a function of observed NEE
plot(O,crps)


# ## CART residual mining
# ## define error metric and dependent variables
# err = (E-O)/beta
# x = inputs$PAR[qaqc]
# colnames(x) = c("temp")
# smp = sample.int(length(err),1000)  ## take a sample of the data since some alg. are slow
# 
# ### Classification tree
# rpb = rpart(err ~ x) ## bias
# plot(rpb)
# text(rpb)
# e2 = err^2
# rpe = rpart(e2 ~ x) ## sq error
# plot(rpe)
# text(rpe)
# 
# ## Random Forest
# rfe = randomForest(x[smp,],abs(err[smp]))
# rfe$importance
# partialPlot(rfe,x[smp,],"PAR")
# partialPlot(rfe,x[smp,],"temp")
# 
# 
# ## raw
# plot(inputs$temp[qaqc],O,pch=".",ylab="NEE")
# points(inputs$temp[qaqc],E,pch=".",col=2)
# 
# ## binned
# nbin = 25
# #PAR = inputs$PAR[qaqc]
# #x = seq(min(PAR),max(PAR),length=nbin)
# Tair = inputs$temp[qaqc]
# xd = seq(min(Tair),max(Tair),length=nbin)
# xmid = xd[-length(xd)] + diff(xd)
# bin = cut(Tair,xd)
# Obar = tapply(O,bin,mean,na.rm=TRUE)
# Ose  = tapply(O,bin,std.error,na.rm=TRUE)
# Ebar = tapply(E,bin,mean,na.rm=TRUE)
# Ese  = tapply(E,bin,std.error,na.rm=TRUE)
# OCI = -cbind(Obar-1.96*Ose,Obar,Obar+1.96*Ose)
# ECI = -cbind(Ebar-1.96*Ese,Ebar,Ebar+1.96*Ese)
# rng = range(rbind(OCI,ECI))
# 
# col2=ecoforecastR::col.alpha("darkgrey",0.9)
# col1=ecoforecastR::col.alpha("lightgrey",0.6)
# 
# plot(xmid,Obar,ylim=rng,type='n',xlab="Air Temperature (C)",ylab="NEP (umol/m2/s)",cex.lab=1.3)
# ecoforecastR::ciEnvelope(xmid,ECI[,1],ECI[,3],col=col2)
# lines(xmid,ECI[,2],col="white",lwd=4)
# ecoforecastR::ciEnvelope(xmid,OCI[,1],OCI[,3],col=col1)
# lines(xmid,OCI[,2],col="lightgrey",lwd=4)
# 
# legend("bottom",legend=c("Model","Data"),lwd=10,col=c(col2,col1),lty=1,cex=1.7)
# 
# 
# ### other summary figures to go in multi-panel
# par(mfrow=c(2,2))
# 
# ## Time-series visualization, daily means
# DoY = floor(L4$DoY-0.02)
# uDoY = sort(unique(DoY))
# ci.pf  = apply(apply(NEE.pf[,],2,tapply,DoY,mean),1,mean)
# NEE = -L4$NEE_st_fMDS
# NEEd = tapply(NEE,DoY,mean)
# plot(uDoY,ci.pf,xlab="time",ylab="NEE",type='l',ylim=range(c(ci.pf,NEEd)),cex.lab=1.3)
# points(uDoY,NEEd,col=2,pch="+")
# legend("topright",legend=c("Model","Data"),lty=c(1,NA),pch=c(NA,"+"),col=1:2,cex=1.3)
# 
# ## predicted vs observed
# plot(NEEd,ci.pf,xlab="Model",ylab="Data",cex.lab=1.3)
# abline(0,1,lty=2,lwd=4)
# abline(lm(ci.pf ~ NEEd),col=2,lwd=3,lty=3)
# legend("topleft",legend=c("1:1","Reg"),lty=2:3,lwd=4,col=1:2,cex=1.3)
# 
# ## Functional response
# plot(xmid,Obar,ylim=rng,type='n',xlab="Air Temperature (C)",ylab="NEP (umol/m2/s)",cex.lab=1.3)
# ecoforecastR::ciEnvelope(xmid,ECI[,1],ECI[,3],col=col2)
# lines(xmid,ECI[,2],col="white",lwd=4)
# ecoforecastR::ciEnvelope(xmid,OCI[,1],OCI[,3],col=col1)
# lines(xmid,OCI[,2],col="lightgrey",lwd=4)
# 
# legend("bottom",legend=c("Model","Data"),lwd=10,col=c(col2,col1),lty=1,cex=1.3)
# 
# ### Classification tree
# par(mar=c(0,0,0,0))
# rpe = rpart(e2 ~ PAR+temp,as.data.frame(x),method="anova") ## sq error
# plot(rpe,margin=0.1)
# text(rpe,cex=1.5)
# 
# 
# 
# ## flux "climatology"
# fluxfiles = dir("data",pattern="AMF")
# fluxfiles = fluxfiles[grep("txt",fluxfiles)]
# fluxfiles = fluxfiles[-grep("2005",fluxfiles)]
# clim.NEE = clim.doy = NULL
# for(f in fluxfiles){
#   ff = read.csv(file.path("data",f),header=TRUE,na.strings="-9999")
#   ff[ff == -9999] = NA
#   clim.NEE = c(clim.NEE,ff$NEE_st_fMDS)
#   clim.doy = c(clim.doy,ff$DoY)
# }
# NEE.clim=tapply(clim.NEE,clim.doy,mean,na.rm=TRUE)[1:length(qaqc)]
# C = NEE.clim[qaqc]
# NEE.clim.fit = lm(O ~ C)
# summary(NEE.clim.fit)
# stats["RMSE",3]  = sqrt(mean((C-O)^2))
# stats['Bias',3]  = mean(C-O)
# stats['cor',3]   = cor(C,O)
# stats['slope',3] = coef(NEE.clim.fit)[2]
# colnames(stats)[3] <- "clim"
# knitr::kable(stats)
# plot(C,O,pch=".",xlab="climatology",ylab='observed',main='NEE (umol/m2/sec)')
# abline(0,1,col=2,lwd=2)
# abline(NEE.clim.fit,col=3,lwd=3,lty=2)
# legend("bottomright",legend=c('obs','1:1','reg'),col=1:3,lwd=3)
# 
# ## example cycle
# plot(L4$DoY,-L4$NEE_st_fMDS,xlim=c(200,210),type='l',lwd=2,ylim=c(-10,20),xlab="Day of Year",ylab="NEE")
# lines(L4$DoY,-NEE.clim,col=4,lwd=2,lty=2)
# legend("topright",legend=c("Obs","clim"),lty=1:2,col=c(1,4),lwd=2)
# 
# 
# ## diurnal cycle
# NEE.ens.diurnal = tapply(E,L4$Hour[qaqc],mean)
# NEE.pf.diurnal  = tapply(P,L4$Hour[qaqc],mean)
# NEE.clim.diurnal  = tapply(C,L4$Hour[qaqc],mean)
# NEE.obs.diurnal = tapply(O,L4$Hour[qaqc],mean)
# ylim=range(c(NEE.ens.diurnal,NEE.pf.diurnal,NEE.obs.diurnal))
# tod = sort(unique(L4$Hour))
# plot(tod,NEE.ens.diurnal,ylim=ylim,col=2,xlab="Time of Day",ylab='NEE',main="Diurnal Cycle",type='l',lwd=3)
# lines(tod,NEE.pf.diurnal,col=3,lwd=3)
# lines(tod,NEE.clim.diurnal,col=4,lwd=3)
# lines(tod,NEE.obs.diurnal,lwd=3)
# legend("bottomright",legend=c("obs","ens","PF","clim"),col=1:4,pch=20,cex=0.75)
