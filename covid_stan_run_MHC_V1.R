###### ######  SETTINGS  ###### ###### ###### 
rm(list=ls())
setwd("~/Desktop/COVID/nowcasting")
mTrsp <- function(cl,a)  { apply(col2rgb(cl), 2, function(x){ rgb(x[1],x[2],x[3],a,maxColorValue=255)}) }

library(Hmisc)
library(splines)
library(rstan)
library(lubridate)
library(RColorBrewer)
library(splines)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

###### ######  INPUTS  ###### ###### ###### 

### Dummy data
set.seed(123)
zz <- round(rgamma(1000,5,0.3))
diagnosis_day0 <- max(zz)-zz+1
# hist(diagnosis_day0)

set.seed(124)
zz <- round(rgamma(1000,3,.9))
days_delay0 <- zz

rmv <- (diagnosis_day0+days_delay0)>max(diagnosis_day0)

days_delay <- days_delay0[!rmv]
diagnosis_day <- diagnosis_day0[!rmv]

days_extra <- 1
# hist(days_delay0); hist(days_delay,add=T,col=5)
# hist(diagnosis_day0); hist(diagnosis_day,add=T,col=5)
range(diagnosis_day)
### Some real data
# raw_data <- read.csv("cases_20200325_NYC.csv")
# 
# diagnosis_day0 <- as.numeric(mdy(as.character(raw_data$diagnosis_date)))
# report_day0    <- as.numeric(mdy(as.character(raw_data$event_create_date)))
# days_delay     <- report_day0 - diagnosis_day0
# 
# diagnosis_day  <- diagnosis_day0 - min(diagnosis_day0) + 1
# 
# # remove the current day, as seems incomplete
# rmv <- max(diagnosis_day) - diagnosis_day - days_delay == 0
# diagnosis_day <- diagnosis_day[!rmv]
# days_delay    <- days_delay[!rmv]
# days_extra <- 1

# Create reporting triangle
N_days_before <- 10
rep_tri_conf_cases <- matrix(0, max(diagnosis_day) + N_days_before,max(days_delay)+1)

for(i in 1:length(diagnosis_day)){
  rep_tri_conf_cases[(diagnosis_day + N_days_before)[i], days_delay[i]+1] = 
    rep_tri_conf_cases[(diagnosis_day + N_days_before)[i], days_delay[i]+1] + 1;
}

n_spl_par <- 10
des_mat <- bs(1:(N_days+days_extra), df=n_spl_par, degree=3, intercept=T)
# this produces a cubic b-spline with n_spl_par basis functions
spl_basis <- as.matrix(as.data.frame(des_mat))

## creating priors from Lauer paper
# progression from infected to symptomatic 
func1 <- function(zz){
  z <- exp(zz)
  prd <- qgamma(c(1,20,39)/40,1/z[1]/z[1],1/z[1]/z[1]/z[2])
  sum( (prd-c(2.2,5.1,11.5))^2 * c(1,5,1))  }

opt_par1 <- exp(optim(c(-1,-1),func1,method="BFGS")$par)
qgamma(c(1,20,39)/40,opt_par1[1]^-2,opt_par1[1]^-2/opt_par1[2])

func2 <- function(z){
  prd <- qgamma(c(1,39)/40,1/z[1]/z[1],1/z[1]/z[1]/5.1)
  sum( (prd-c(4.5,5.8))^2)   }

opt_par2 <- optimize(func2,c(0.01,1))$minimum
qgamma(c(1,39)/40,opt_par2^-2,opt_par2^-2/5.1)

###### ######  LIST  ###### ###### ###### 

datList <- list()

## Data
datList[["N_conf_cases"]]     <-  length(diagnosis_day) #
datList[["N_days"]]           <-  max(diagnosis_day) + N_days_before #n days to model before first case
datList[["N_days_extra"]]     <-  days_extra # model days after (for things like death, hospitalizations)
datList[["Max_delay"]]        <-  max(days_delay) # maximum value for delay distribution
datList[["cases_test_day"]]   <-  diagnosis_day + N_days_before
datList[["cases_days_delay"]] <-  days_delay # NEED A PRIOR ON THIS

## Priors
# priors on spline
datList[["spl_basis"]]  <-  spl_basis
datList[["n_spl_par"]]  <-  n_spl_par
# probability of transitioning between states 
# unlike ODE, not a constant rate in transmission -- allows us to model delays
# probability of transitioning into the next stage or recovery
datList[["pri_p_sym_if_inf_mn"]] <-  0.69
datList[["pri_p_sym_if_inf_ss"]] <-  100 # sum of the beta coeffecients

datList[["pri_p_hos_if_sym_mn"]] <-  0.31 # upper bound from CDC paper
datList[["pri_p_hos_if_sym_ss"]] <-  100

datList[["pri_p_die_if_hos_mn"]] <-  0.03 # upper bound from CDC paper
datList[["pri_p_die_if_hos_ss"]] <-  100

datList[["nb_yes"]]              <-  0

#### Delay from infection to symptoms
# https://annals.org/aim/fullarticle/2762808/incubation-period-coronavirus-disease-2019-covid-19-from-publicly-reported
#  Median incubation 5.1 days (4.5, 5.8 ); 97.5 percentile 11.5 days (8.2, 15.6); 2.5 percentile 2.2 days (1.8, 2.9) 
# functions to get parameters from the Lauer paper

# shape is fixed, mean can vary
datList[["pri_inf_prg_delay_mn_mn"]] <-  opt_par1[2] # Lauer et al. 
datList[["pri_inf_prg_delay_mn_cv"]] <-  opt_par2
datList[["inf_prg_delay_cv"]]        <-  opt_par1[1]
# prior on the mean of the gamma distribution is distributed gamma with a mean and a cv
# we could simplify by fixing the mean OR assuming time to recovery and time to progression is same
datList[["pri_sym_prg_delay_mn"]]  <-  11 # Zhou et al. 
datList[["pri_sym_prg_delay_low"]] <-  8
datList[["pri_sym_prg_delay_up"]]      <-  14

datList[["pri_hos_prg_delay_mn"]] <-  8.8 # Linton et al. 
datList[["pri_hos_prg_delay_low"]] <-  7.2
datList[["pri_hos_prg_delay_up"]]  <-  10.8

datList[["pri_inf_res_delay_mn_mn"]] <-  14 # assumption 
datList[["pri_inf_res_delay_mn_cv"]] <-  0.1
datList[["inf_res_delay_cv"]]        <-  0.5

datList[["pri_sym_res_delay_mn_mn"]] <-  6 #assumption
datList[["pri_sym_res_delay_mn_cv"]] <-  0.1
datList[["sym_res_delay_cv"]]        <-  0.5

datList[["pri_hos_res_delay_mn"]] <-  17 # Pan et al. 
datList[["pri_hos_res_delay_low"]] <-  13
datList[["pri_hos_res_delay_up"]]  <-  21

datList[["pri_report_delay_mn_mn"]]  <-  7
datList[["pri_report_delay_mn_cv"]]  <-  0.9

datList[["pri_report_delay_cv_mn"]]  <-  0.5
datList[["pri_report_delay_cv_cv"]]  <-  0.9
# daily probability of diagnosis as a function of individual health state
# currently parameterizing beta distributions, but we need to add a function of time
    # probability of diagnosis is a logistic function with intercepts by health state
datList[["pri_p_diag_if_inf_mn"]] <-  0.01
datList[["pri_p_diag_if_inf_ss"]] <-  100
datList[["pri_p_diag_if_sym_mn"]] <-  0.1
datList[["pri_p_diag_if_sym_ss"]] <-  100
datList[["pri_p_diag_if_hos_mn"]] <-  0.6
datList[["pri_p_diag_if_hos_ss"]] <-  100

###### ######  RUN MODEL  ###### ###### ###### 

source("covid_stan_script_MHC_V1.stan")

fit_stan <-  stan(model_code = stan_code,
                  control=list(adapt_delta=0.92, 
                               max_treedepth=12),
                  data = datList, 
                  seed = 1234, 
                  chains = 3, 
                  iter = 500, 
                  warmup = 400)

samps <- extract(fit_stan)

# summary(fit_stan)
# names(samps)

###### ######  PLOTS  ###### ###### ###### 
smp  <- ceiling(seq(1,nrow( samps[["new_inf"]]),length.out=10))

########### Comparison to data ########### 
pdfnam <- "Fig_COVID_nowcast_temp_4-7-2020.pdf"
pdf(file=pdfnam,width=9, height=6.5) 

N_days <- datList[["N_days"]]
N_days_tot <- datList[["N_days"]] + datList[["N_days_extra"]]

max_delay <- max(days_delay)+1
cls <- colorRampPalette(c("red4","red","blue","blue4"))(max_delay)

otcm_mu <- apply(samps[["rep_tri_conf_cases_mu"]],2:3,mean)

par(mfrow=c(4,3), mar=c(1.2,1.5,.2,.2),oma=c(0.5,3,3.5,.5))
for(i in 1:(max_delay)){
  plot(1:(N_days-i+1),otcm_mu[1:(N_days-i+1),i],ylim=c(0,max(c(otcm_mu[1:(N_days-i+1),i],rep_tri_conf_cases[1:(N_days-i+1),i])))*1.12,xlim=c(N_days_before,N_days),type="l",col=cls[i],
       axes=F,xlab="",ylab=""); 
  axis(2,las=1,tcl=-.07,mgp=c(3, 0.15, 0));axis(1,las=1,tcl=-.07,mgp=c(3, 0.1, 0));box()
  text(N_days_before,max(c(otcm_mu[1:(N_days-i+1),i],rep_tri_conf_cases[1:(N_days-i+1),i]))*1.06,paste0(i-1," days delay"),pos=4,offset=-.2,font=2,cex=1.15)
  lines(1:(N_days-i+1),rep_tri_conf_cases[1:(N_days-i+1),i],col=cls[i],lty="12");points(1:(N_days-i+1),rep_tri_conf_cases[1:(N_days-i+1),i],col=cls[i],pch=16,cex=.6)
}
mtext("Case count (N)",2,1,outer=T)
mtext("Fit to data: Case time-series, stratified by days delay (points = data, lines = model)",3,1,outer=T)

########### Comparison to data 2 ########### 

par(mfrow=c(1,1), mar=c(2.5,3.5,2.0,.2),oma=c(0,0,0,0))

otcm_rep_del <- samps[["rep_tri_conf_cases_mu"]]
for(k in 1:dim(otcm_rep_del)[1]) { for(i in 1:N_days) {  for(j in 1:max_delay) {
  if((i+j) >= (N_days+2)){ otcm_rep_del[k,i,j] = 0 ;  }
}  }  }

otcm_interval <- apply(apply(samps[["rep_tri_conf_cases_mu"]],1:2,sum),2,function(x) quantile(x,c(1,39)/40))
otcm_interval2 <- apply(apply(otcm_rep_del,1:2,sum),2,function(x) quantile(x,c(1,39)/40))
otcm_mu <- apply(samps[["rep_tri_conf_cases_mu"]],2:3,mean)
otcm_mu2 <- apply(otcm_rep_del,2:3,mean)

plot(1:N_days,rowSums(otcm_mu),axes=F,xlab="",ylab="",type="l",lty="11",col=NA,ylim=c(0,max(otcm_interval)),xlim=c(which(rowSums(rep_tri_conf_cases)>0)[1],N_days))
polygon(c(1:N_days,N_days:1),c(otcm_interval[1,],otcm_interval[2,N_days:1]),border=mTrsp("red",30),col=mTrsp("red",60))
polygon(c(1:N_days,N_days:1),c(otcm_interval2[1,],otcm_interval2[2,N_days:1]),border=mTrsp("navy",30),col=mTrsp("navy",60))
lines(1:N_days,rowSums(otcm_mu),col="red")
lines(1:N_days,rowSums(otcm_mu2),col="navy")
lines(1:N_days,rowSums(rep_tri_conf_cases),col=mTrsp("navy",150),lty="12");
points(1:N_days,rowSums(rep_tri_conf_cases),col="navy",pch=16,cex=.45)
axis(2,las=1,tcl=-.07,mgp=c(3, 0.15, 0));axis(1,las=1,tcl=-.07,mgp=c(3, 0.1, 0));box()
mtext("Case count (N)",2,2.3)
mtext("Days since start",1,1.3)
mtext("Fit to data: Cases diagnosed per day",3,0.6)
legend("topleft",c("Model: new diagnoses (all)","Model: new diagnoses (reported)","Data: new diagnoses (reported)"),
       col=c("red","navy","navy"),lty=c(1,1,3),
       pch=c(NA,NA,16),pt.cex=0.8,seg.len=1.5,lwd=c(2,1,1),cex=0.9,bty="n",ncol=1)
#####
new_inf  <- apply(samps[["new_inf"]],2,mean)
new_sym  <- apply(samps[["new_sym"]],2,mean)
new_hos  <- apply(samps[["new_hos"]],2,mean)
new_die  <- apply(samps[["new_die"]],2,mean)
new_res  <- apply(samps[["res_inf"]]+samps[["res_sym"]]+samps[["res_hos"]],2,mean)
mx <- max(new_inf,new_sym,new_die,new_res)

plot(1:N_days,rowSums(otcm_mu),axes=F,xlab="",ylab="",type="l",lty="11",col=NA,ylim=c(1,mx*1.1),xlim=c(which(rowSums(rep_tri_conf_cases)>0)[1],N_days_tot),log="y")
lines(1:N_days_tot,new_inf,col="grey40",lwd=2)
lines(1:N_days_tot,new_sym,col=2,lwd=2)
lines(1:N_days_tot,new_hos,col="forestgreen",lwd=2)
lines(1:N_days_tot,new_die,col="purple",lwd=2)
lines(1:N_days_tot,new_res,col="brown",lwd=2)
lines(1:N_days,rowSums(otcm_mu),col="orange",lwd=2)
lines(1:N_days,rowSums(otcm_mu2),col="navy",lwd=2)
lines(1:N_days,rowSums(rep_tri_conf_cases),col=mTrsp("navy",150),lty="12");
points(1:N_days,rowSums(rep_tri_conf_cases),col="navy",pch=16,cex=.6)
axis(2,las=1,tcl=-.07,mgp=c(3, 0.15, 0));axis(1,las=1,tcl=-.07,mgp=c(3, 0.1, 0));box()
mtext("Case count (N)",2,2.3)
mtext("Days since start",1,1.3)
mtext("Fit to data: Incident outcomes compared to data",3,0.6)
legend("topleft",c("Model: new infections","Model: new symptomatics","Model: new hospitalizations",
       "Model: deaths","Model: recoveries",
       "Model: new diagnoses (all)","Model: new diagnoses (reported)","Data: new diagnoses (reported)"),
       col=c("grey40","red","forestgreen","blue","purple","brown","orange","navy","navy"),lty=c(rep(1,8),NA),
       pch=c(rep(NA,8),16),pt.cex=0.7,seg.len=1.2,lwd=2,cex=0.9,bty="n",ncol=2)

########### Incident Outcomes ########### 
#### New Infections
par(mfrow=c(2,3), mar=c(1.2,2.0,.2,.2),oma=c(3.5,3,3.5,.5))

otcm  <- samps[["new_inf"]]
plot(1:N_days_tot,apply(otcm,2,mean),type="l",ylim=c(0,quantile(otcm[,N_days_tot],.975))*1.04,las=1,
     main="",axes=F,xlab="",ylab="")
axis(2,las=1,tcl=-.07,mgp=c(3, 0.15, 0));axis(1,las=1,tcl=-.07,mgp=c(3, 0.1, 0));box()
polygon(c(1:N_days_tot,N_days_tot:1),c(apply(otcm,2,function(x) quantile(x,0.025)),
                                 apply(otcm,2,function(x) quantile(x,0.975))[N_days_tot:1]),border="grey82",col="grey90")
for(i in smp) lines(1:N_days_tot,otcm[i,],col=mTrsp(3,120))
lines(1:N_days_tot,apply(otcm,2,mean))
text(0,quantile(otcm[,N_days_tot],.975)*1.01,"New infections per day",pos=4,offset=0.2,font=2,cex=1.15)

#### New Symp
otcm <- samps[["new_sym"]]
plot(1:N_days_tot,apply(otcm,2,mean),type="l",ylim=c(0,quantile(otcm[,N_days_tot],.975))*1.04,las=1,
     main="",xlab="",ylab="",axes=F)
axis(2,las=1,tcl=-.07,mgp=c(3, 0.15, 0));axis(1,las=1,tcl=-.07,mgp=c(3, 0.1, 0));box()

polygon(c(1:N_days_tot,N_days_tot:1),c(apply(otcm,2,function(x) quantile(x,0.025)),
                             apply(otcm,2,function(x) quantile(x,0.975))[N_days_tot:1]),border="grey82",col="grey90")
for(i in smp) lines(1:N_days_tot,otcm[i,],col=mTrsp(3,120))
lines(1:N_days_tot,apply(otcm,2,mean))
text(0,quantile(otcm[,N_days_tot],.975)*1.01,"New symptomatic per day",pos=4,offset=.2,font=2,cex=1.15)

#### New Hosp
otcm <- samps[["new_hos"]]
plot(1:N_days_tot,apply(otcm,2,mean),type="l",ylim=c(0,quantile(otcm[,N_days_tot],.975))*1.04,las=1,
     main="",xlab="",ylab="",axes=F)
axis(2,las=1,tcl=-.07,mgp=c(3, 0.15, 0));axis(1,las=1,tcl=-.07,mgp=c(3, 0.1, 0));box()

polygon(c(1:N_days_tot,N_days_tot:1),c(apply(otcm,2,function(x) quantile(x,0.025)),
                             apply(otcm,2,function(x) quantile(x,0.975))[N_days_tot:1]),border="grey82",col="grey90")
for(i in smp) lines(1:N_days_tot,otcm[i,],col=mTrsp(3,120))
lines(1:N_days_tot,apply(otcm,2,mean))
text(0,quantile(otcm[,N_days_tot],.975)*1.01,"New hospitalizations per day",pos=4,offset=.2,font=2,cex=1.15)

#### New death
otcm <- samps[["new_die"]]
plot(1:N_days_tot,apply(otcm,2,mean),type="l",ylim=c(0,quantile(otcm[,N_days_tot],.975))*1.04,las=1,
     main="",xlab="",ylab="",axes=F)
axis(2,las=1,tcl=-.07,mgp=c(3, 0.15, 0));axis(1,las=1,tcl=-.07,mgp=c(3, 0.1, 0));box()

polygon(c(1:N_days_tot,N_days_tot:1),c(apply(otcm,2,function(x) quantile(x,0.025)),
                             apply(otcm,2,function(x) quantile(x,0.975))[N_days_tot:1]),border="grey82",col="grey90")
for(i in smp) lines(1:N_days_tot,otcm[i,],col=mTrsp(3,120))
lines(1:N_days_tot,apply(otcm,2,mean))
text(0,quantile(otcm[,N_days_tot],.975)*1.01,"Deaths per day",pos=4,offset=.2,font=2,cex=1.15)

#### New recovery
otcm <- samps[["res_inf"]]+samps[["res_sym"]]+samps[["res_hos"]]
plot(1:N_days_tot,apply(otcm,2,mean),type="l",ylim=c(0,quantile(otcm[,N_days_tot],.975))*1.04,las=1,
     main="",xlab="",ylab="",axes=F)
axis(2,las=1,tcl=-.07,mgp=c(3, 0.15, 0));axis(1,las=1,tcl=-.07,mgp=c(3, 0.1, 0));box()

polygon(c(1:N_days_tot,N_days_tot:1),c(apply(otcm,2,function(x) quantile(x,0.025)),
                             apply(otcm,2,function(x) quantile(x,0.975))[N_days_tot:1]),border="grey82",col="grey90")
for(i in smp) lines(1:N_days_tot,otcm[i,],col=mTrsp(3,120))
lines(1:N_days_tot,apply(otcm,2,mean))
text(0,quantile(otcm[,N_days_tot],.975)*1.01,"Recoveries per day",pos=4,offset=.2,font=2,cex=1.15)

mtext("Count of outcomes (N)",2,1,outer=T)
mtext("Incident outcomes: new infections, symptomatics, hospitalizations, deaths, and recoveries",3,1,outer=T)
mtext("Days since start",1,1,outer=T)


########### Current outcomes and cumulative final outcomes ########### 
par(mfrow=c(2,3), mar=c(1.2,2.1,.2,.2),oma=c(3.5,3,3.5,.5))

#### Current Asymptomatic
otcm <- samps[["cur_inf"]]
plot(1:N_days_tot,apply(otcm,2,mean),type="l",ylim=c(0,quantile(otcm[,N_days_tot],.975))*1.04,las=1,
     main="",xlab="",ylab="",axes=F)
axis(2,las=1,tcl=-.07,mgp=c(3, 0.15, 0));axis(1,las=1,tcl=-.07,mgp=c(3, 0.1, 0));box()

polygon(c(1:N_days_tot,N_days_tot:1),c(apply(otcm,2,function(x) quantile(x,0.025)),
                             apply(otcm,2,function(x) quantile(x,0.975))[N_days_tot:1]),border="grey82",col="grey90")
for(i in smp) lines(1:N_days_tot,otcm[i,],col=mTrsp(3,120))
lines(1:N_days_tot,apply(otcm,2,mean))
text(0,quantile(otcm[,N_days_tot],.975)*1.01,"Current asymptomatic",pos=4,offset=.2,font=2,cex=1.15)

#### Current Symptomatic, Non-Hospitalized
otcm <- samps[["cur_sym"]]
plot(1:N_days_tot,apply(otcm,2,mean),type="l",ylim=c(0,quantile(otcm[,N_days_tot],.975))*1.04,las=1,
     main="",xlab="",ylab="",axes=F)
axis(2,las=1,tcl=-.07,mgp=c(3, 0.15, 0));axis(1,las=1,tcl=-.07,mgp=c(3, 0.1, 0));box()

polygon(c(1:N_days_tot,N_days_tot:1),c(apply(otcm,2,function(x) quantile(x,0.025)),
                             apply(otcm,2,function(x) quantile(x,0.975))[N_days_tot:1]),border="grey82",col="grey90")
for(i in smp) lines(1:N_days_tot,otcm[i,],col=mTrsp(3,120))
lines(1:N_days_tot,apply(otcm,2,mean))
text(0,quantile(otcm[,N_days_tot],.975)*1.01,"Current symptomatic, non-hospitalized",pos=4,offset=.2,font=2,cex=1.15)

#### Current Hospitalized
otcm <- samps[["cur_hos"]]
plot(1:N_days_tot,apply(otcm,2,mean),type="l",ylim=c(0,quantile(otcm[,N_days_tot],.975))*1.04,las=1,
     main="",xlab="",ylab="",axes=F)
axis(2,las=1,tcl=-.07,mgp=c(3, 0.15, 0));axis(1,las=1,tcl=-.07,mgp=c(3, 0.1, 0));box()

polygon(c(1:N_days_tot,N_days_tot:1),c(apply(otcm,2,function(x) quantile(x,0.025)),
                             apply(otcm,2,function(x) quantile(x,0.975))[N_days_tot:1]),border="grey82",col="grey90")
for(i in smp) lines(1:N_days_tot,otcm[i,],col=mTrsp(3,120))
lines(1:N_days_tot,apply(otcm,2,mean))
text(0,quantile(otcm[,N_days_tot],.975)*1.01,"Current hospitalized",pos=4,offset=.2,font=2,cex=1.15)

#### Cumulative dead
otcm <- samps[["cum_die"]]
plot(1:N_days_tot,apply(otcm,2,mean),type="l",ylim=c(0,quantile(otcm[,N_days_tot],.975))*1.04,las=1,
     main="",xlab="",ylab="",axes=F)
axis(2,las=1,tcl=-.07,mgp=c(3, 0.15, 0));axis(1,las=1,tcl=-.07,mgp=c(3, 0.1, 0));box()

polygon(c(1:N_days_tot,N_days_tot:1),c(apply(otcm,2,function(x) quantile(x,0.025)),
                             apply(otcm,2,function(x) quantile(x,0.975))[N_days_tot:1]),border="grey82",col="grey90")
for(i in smp) lines(1:N_days_tot,otcm[i,],col=mTrsp(3,120))
lines(1:N_days_tot,apply(otcm,2,mean))
text(0,quantile(otcm[,N_days_tot],.975)*1.01,"Cumulative dead",pos=4,offset=.2,font=2,cex=1.15)

#### Cumulative recovered
otcm <- samps[["cum_res"]]
plot(1:N_days_tot,apply(otcm,2,mean),type="l",ylim=c(0,quantile(otcm[,N_days_tot],.975))*1.04,las=1,
     main="",xlab="",ylab="",axes=F)
axis(2,las=1,tcl=-.07,mgp=c(3, 0.15, 0));axis(1,las=1,tcl=-.07,mgp=c(3, 0.1, 0));box()

polygon(c(1:N_days_tot,N_days_tot:1),c(apply(otcm,2,function(x) quantile(x,0.025)),
                             apply(otcm,2,function(x) quantile(x,0.975))[N_days_tot:1]),border="grey82",col="grey90")
for(i in smp) lines(1:N_days_tot,otcm[i,],col=mTrsp(3,120))
lines(1:N_days_tot,apply(otcm,2,mean))
text(0,quantile(otcm[,N_days_tot],.975)*1.01,"Cumulative recovered",pos=4,offset=.2,font=2,cex=1.15)

mtext("Count of outcomes (N)",2,1,outer=T)
mtext("Current and cumulative outcomes: new infections, symptomatics, hospitalizations, deaths, and recoveries",3,1,outer=T)
mtext("Days since start",1,1,outer=T)

############ Priors vs. Posteriors ########### 

par(mfrow=c(3,4), mar=c(1.2,1.3,2,.2),oma=c(0.5,2.5,.5,.5))

##  log_new_inf_0
mu = datList[["pri_log_new_inf_0_mu"]]; sd = datList[["pri_log_new_inf_0_sd"]]
hist(samps[["log_new_inf_0"]],col="goldenrod",border=F,breaks=20,axes=F,main="",xlab="",ylab="",xlim=range(c(samps[["log_new_inf_0"]],mu-sd/2,mu+sd/2)),probability = T)
axis(2,las=1,tcl=-.07,mgp=c(3, 0.15, 0));axis(1,las=1,tcl=-.07,mgp=c(3, 0.1, 0));box()
zz <- range(c(samps[["log_new_inf_0"]],mu-sd,mu+sd))
x_pts <- seq(zz[1],zz[2],length.out=101)
y_pts <- dnorm(x_pts,mu,sd)
lines(x_pts,y_pts,col=mTrsp(4,200),lwd=2)
mtext("log_new_inf_0",3,0.5,font=2,cex=0.8)
##  log_new_inf_drift
mu = datList[["pri_log_new_inf_drift_mu"]]; sd = datList[["pri_log_new_inf_drift_sd"]]
mod = samps[["log_new_inf_drift"]]
hist(mod,col="goldenrod",border=F,breaks=10,axes=F,main="",xlab="",ylab="",xlim=range(c(mod,mu-sd/2,mu+sd/2)),probability = T)
axis(2,las=1,tcl=-.07,mgp=c(3, 0.15, 0));axis(1,las=1,tcl=-.07,mgp=c(3, 0.1, 0));box()
zz <- range(c(mod,mu-sd,mu+sd))
x_pts <- seq(zz[1],zz[2],length.out=101)
y_pts <- dnorm(x_pts,mu,sd)
lines(x_pts,y_pts,col=mTrsp(4,200),lwd=2)
mtext("log_new_inf_drift",3,0.5,font=2,cex=0.8)
##  pri_sigma_deriv1_log_new_inf_sd
mu = 0; sd = datList[["pri_sigma_deriv1_log_new_inf_sd"]]
mod = samps[["sigma_deriv1_log_new_inf"]]
hist(mod,col="goldenrod",border=F,breaks=20,axes=F,main="",xlab="",ylab="",xlim=range(c(mod,0,mu+sd/2)),probability = T)
axis(2,las=1,tcl=-.07,mgp=c(3, 0.15, 0));axis(1,las=1,tcl=-.07,mgp=c(3, 0.1, 0));box()
zz <- range(c(mod,mu-sd,mu+sd))
x_pts <- seq(zz[1],zz[2],length.out=101)
y_pts <- dnorm(x_pts,mu,sd)
lines(x_pts,y_pts,col=mTrsp(4,200),lwd=2)
mtext("sigma_deriv1_log_new_inf",3,0.5,font=2,cex=0.8)
## p_sym_if_inf
mn = datList[["pri_p_sym_if_inf_mn"]]; ss = datList[["pri_p_sym_if_inf_ss"]]
mod = samps[["p_sym_if_inf"]]
hist(mod,col="goldenrod",border=F,breaks=15,axes=F,main="",xlab="",ylab="",xlim=range(c(mod,qbeta(c(0.0002,0.9998),mn*ss,(1-mn)*ss))),probability = T)
axis(2,las=1,tcl=-.07,mgp=c(3, 0.15, 0));axis(1,las=1,tcl=-.07,mgp=c(3, 0.1, 0));box()
zz <- range(c(mod,0,1))
x_pts <- seq(zz[1],zz[2],length.out=101)
y_pts <- dbeta(x_pts,mn*ss,(1-mn)*ss)
lines(x_pts,y_pts,col=mTrsp(4,200),lwd=2)
mtext("p_sym_if_inf",3,0.5,font=2,cex=0.8)
## p_hos_if_sym
mn = datList[["pri_p_hos_if_sym_mn"]]; ss = datList[["pri_p_hos_if_sym_ss"]]
mod = samps[["p_hos_if_sym"]]
hist(mod,col="goldenrod",border=F,breaks=15,axes=F,main="",xlab="",ylab="",xlim=range(c(mod,qbeta(c(0.0002,0.9998),mn*ss,(1-mn)*ss))),probability = T)
axis(2,las=1,tcl=-.07,mgp=c(3, 0.15, 0));axis(1,las=1,tcl=-.07,mgp=c(3, 0.1, 0));box()
zz <- range(c(mod,0,1))
x_pts <- seq(zz[1],zz[2],length.out=101)
y_pts <- dbeta(x_pts,mn*ss,(1-mn)*ss)
lines(x_pts,y_pts,col=mTrsp(4,200),lwd=2)
mtext("p_hos_if_sym",3,0.5,font=2,cex=0.8)
## p_die_if_hos
mn = datList[["pri_p_die_if_hos_mn"]]; ss = datList[["pri_p_die_if_hos_ss"]]
mod = samps[["p_die_if_hos"]]
hist(mod,col="goldenrod",border=F,breaks=15,axes=F,main="",xlab="",ylab="",xlim=range(c(mod,qbeta(c(0.0002,0.9998),mn*ss,(1-mn)*ss))),probability = T)
axis(2,las=1,tcl=-.07,mgp=c(3, 0.15, 0));axis(1,las=1,tcl=-.07,mgp=c(3, 0.1, 0));box()
zz <- range(c(mod,0,1))
x_pts <- seq(zz[1],zz[2],length.out=101)
y_pts <- dbeta(x_pts,mn*ss,(1-mn)*ss)
lines(x_pts,y_pts,col=mTrsp(4,200),lwd=2)
mtext("p_die_if_hos",3,0.5,font=2,cex=0.8)

mtext("Model parameters: prior (blue) versus posterior (yellow)",3,1,outer=T)

dev.off();system(paste("open", pdfnam))


names(samps)
names(datList)

########################


