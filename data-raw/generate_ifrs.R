#!/usr/bin/Rscript

ifrs <- local({
  # ~~~~~~~~~~~~~~~~~~~~~~ [1] obtain and clean death data ~~~~~~~~~~~~~~~~~~~~~~

  # get data from https://data.cdc.gov/NCHS/Provisional-COVID-19-Death-Counts-by-Sex-Age-and-S/9bhg-hcku
  age_deaths <- read.csv("data-raw/ifr-data/Provisional_COVID-19_Death_Counts_by_Sex__Age__and_State.csv")

  # remove overlapping age-groups
  kp_ag1 <- unique(age_deaths$Age.group)[c(2,4,5,6,8,10,11,13,14,15,16)]
  kp_ag2 <- age_deaths$Age.group%in%kp_ag1

  # remove 'all' sex group, as that only available at nat level.
  kp_sx1 <- unique(age_deaths$Sex)[2:3]
  kp_sx2 <- age_deaths$Sex%in%kp_sx1
  age_deaths1 <-  age_deaths[kp_ag2 & kp_sx2,4:7]

  # implement a simple mean imputation for interval-censored data
  tot_death1 <- sum(age_deaths1$COVID.19.Deaths[age_deaths1$State=="United States"])
  tot_death2 <- sum(age_deaths1$COVID.19.Deaths[age_deaths1$State!="United States"],na.rm=T)
  imp_death_n <- (tot_death1-tot_death2)/sum(is.na(age_deaths1$COVID.19.Deaths))
  age_deaths1$COVID.19.Deaths_imp <- age_deaths1$COVID.19.Deaths
  age_deaths1$COVID.19.Deaths_imp[is.na(age_deaths1$COVID.19.Deaths)] <- imp_death_n

  # remove puerto rico and United States total, and add NYC into NY state
  age_deaths1$COVID.19.Deaths_imp[age_deaths1$State=="New York"] <- age_deaths1$COVID.19.Deaths_imp[age_deaths1$State=="New York"]+
                                                                    age_deaths1$COVID.19.Deaths_imp[age_deaths1$State=="New York City"]
  age_deaths1 <- age_deaths1[!age_deaths1$State%in%c("New York City","Puerto Rico","United States"),]

  age_deaths2 <- aggregate(COVID.19.Deaths_imp ~State+Age.group,age_deaths1,sum)

  # ~~~~~~~~~~~~~~~~~~~~~~ [2] stan model to smooth death counts by age and state ~~~~~~~~~~~~~~~~~~~~~~

  # Poisson model with state and age-group fixed effects, and random effects for interaction of state and age.
  library(rstan)
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())

  age_deaths2$age_state <- paste0(age_deaths2$State,"_",age_deaths2$Age.group)

  datList <- list(N_obs       = nrow(age_deaths2),
                  N_state     = length(unique(age_deaths2$State)), 
                  N_age       = length(unique(age_deaths2$Age.group)), 
                  obs_die     = round(age_deaths2$COVID.19.Deaths_imp),
                  state       = as.numeric(as.factor(age_deaths2$State)),
                  age         = as.numeric(as.factor(age_deaths2$Age.group)),
                  age_state   = as.numeric(as.factor(age_deaths2$age_state)),
                  sigma_re_sd = 0.5)

  source("data-raw/smooth_mort_state_V1.stan")

  fit <- rstan::stan(
    model_code = stan_code, 
    control = list(adapt_delta = 0.96, max_treedepth = 12),
    data    = datList,
    seed    = 123,
    chains  = 4,
    iter    = 2250,
    warmup  = 2000)

  extract_fit <- extract(fit)

  # save(extract_fit,file ="fit_state_age_model_10-7-2020.rData")
  # load("fit_state_age_model_10-7-2020.rData") # extract_fit

  age_deaths2$COVID.19.Deaths_smooth <- colMeans(extract_fit$pred)

  # ~~~~~~~~~~~~~~~~~~~~~~ [3] create state-specific IFRs ~~~~~~~~~~~~~~~~~~~~~~

  # CDC IFRs by age from https://www.cdc.gov/coronavirus/2019-ncov/hcp/planning-scenarios.html
  ifr_cdc <- c(0.00003,0.0002,0.005,0.054)
  names(ifr_cdc) <-  c("0_19_yrs","20_49_yrs","50_69_yrs","70_plus_yrs")

  # Apply IFRs, with simple average where groups dont overlap
  age_deaths2$ifr_cdc[age_deaths2$Age.group%in%c("Under 1 year","1-4 years","5-14 years")] <- ifr_cdc[1]
  age_deaths2$ifr_cdc[age_deaths2$Age.group=="15-24 years"] <- mean(ifr_cdc[1:2])
  age_deaths2$ifr_cdc[age_deaths2$Age.group%in%c("25-34 years","35-44 years")] <- ifr_cdc[2]
  age_deaths2$ifr_cdc[age_deaths2$Age.group=="45-54 years"] <- mean(ifr_cdc[2:3])
  age_deaths2$ifr_cdc[age_deaths2$Age.group=="55-64 years"] <- ifr_cdc[3]
  age_deaths2$ifr_cdc[age_deaths2$Age.group=="65-74 years"] <- mean(ifr_cdc[3:4])
  age_deaths2$ifr_cdc[age_deaths2$Age.group=="75-84 years"] <- ifr_cdc[4]
  age_deaths2$ifr_cdc[age_deaths2$Age.group=="85 years and over"] <- ifr_cdc[4]

  age_deaths2$implied_infections <- age_deaths2$COVID.19.Deaths_smooth/age_deaths2$ifr_cdc

  tot_inf <- sum(age_deaths2$implied_infections)/1e6
  tot_death <- sum(age_deaths2$COVID.19.Deaths_smooth)/1e6

  ifr_state <- aggregate(implied_infections~State,sum,data=age_deaths2)
  ifr_state$deaths <- aggregate(COVID.19.Deaths_smooth~State,sum,data=age_deaths2)[,2]
  ifr_state$ifr <- ifr_state$deaths/ifr_state$implied_infections

  # ~~~~~~~~~~~~~~~~~~~~~~ [3] create county-specific IFRs ~~~~~~~~~~~~~~~~~~~~~~

  # load and clean data from https://www.cdc.gov/mmwr/volumes/69/wr/mm6929a1.htm
  ifr_county <- read.csv("data-raw/ifr-data/cdc_90519_DS1.csv")[,c(1,2,5,6,7,10)]
  ifr_county$anycondition_prevalence2 <-  ifr_county$anycondition_number/ifr_county$county_pop2018_18.and.older

  # make state averages, and calculate county average relative to the state
  comorb_state <- aggregate(county_pop2018_18.and.older~STATE_NAME,ifr_county,sum)
  comorb_state$comorb_n <- aggregate(anycondition_number~STATE_NAME,ifr_county,sum)[,2]
  comorb_state$comorb_fract <- comorb_state$comorb_n/comorb_state$county_pop2018_18.and.older

  ifr_county$state_av <- NA
  for(i in 1:nrow(comorb_state)){
    ifr_county$state_av[ifr_county$STATE_NAME==comorb_state$STATE_NAME[i]] <- comorb_state$comorb_fract[i]
  }

  ifr_county$comorb_rel_to_state <- ifr_county$anycondition_prevalence2/ifr_county$state_av

  # calculate county ifrs
  ifr_county$ifr_cty <- NA
  for(i in 1:nrow(ifr_state)){
    ifr_county$ifr_state[ifr_county$STATE_NAME==ifr_state$State[i]] <- ifr_state$ifr[i]
  }
  ifr_county$ifr_cty <- ifr_county$ifr_state*ifr_county$comorb_rel_to_state

  # ~~~~~~~~~~~~~~~~~~~~~~ [4] calculate beta parameters from mean and coef of variation ~~~~~~~~~~~~~~~~~~~~~~

  # Function adapted from here (and checked) https://stats.stackexchange.com/questions/12232/calculating-the-parameters-of-a-beta-distribution-using-the-mean-and-variance
  estBetaPars <- function(mu, cv) {
    alpha <- ((1 - mu) / (cv*mu)^2 - 1 / mu) * mu ^ 2
    beta <- alpha * (1 / mu - 1)
    return( c(alpha,beta) )
  }

  # Implement for CV = 0.25
  # county level
  ifr_county$ifr_par2 <- ifr_county$ifr_par1 <- NA
  for(i in 1: nrow(ifr_county)){
      zz <- estBetaPars(ifr_county$ifr_cty[i],.25)
      ifr_county$ifr_par1[i] <- zz[1]
      ifr_county$ifr_par2[i] <- zz[2]
  }

  # state level
  ifr_state$ifr_par2 <- ifr_state$ifr_par1 <- NA
  for(i in 1: nrow(ifr_state)){
    zz <- estBetaPars(ifr_state$ifr[i],.25)
    ifr_state$ifr_par1[i] <- zz[1]
    ifr_state$ifr_par2[i] <- zz[2]
  }

  # zeros at front of county FIPS
  ifr_county$FIPS <- as.character(ifr_county$FIPS)
  ifr_county$FIPS[nchar(ifr_county$FIPS)==4] <- paste0(0,ifr_county$FIPS[nchar(ifr_county$FIPS)==4])
  ifr_county$FIPS <- as.character(ifr_county$FIPS)
  ifr_county$FIPS[nchar(ifr_county$FIPS)==4] <- paste0(0,ifr_county$FIPS[nchar(ifr_county$FIPS)==4])

  # clean things up a bit
  ifr_state  <- ifr_state[,c(1,4:6)]
  ifr_county <- ifr_county[,c(2,1,3,9,10,12,13)]
  names(ifr_county) <- c("state","county","fips","comorb_rel_to_state","ifr","ifr_par1","ifr_par2")

  list(
    ifr_state = ifr_state,
    ifr_county = ifr_county
  )
});
