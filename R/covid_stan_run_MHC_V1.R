runModel <- function(data = defaultData(), chains=3, iter=500) {
  ###### ###### SETTINGS ###### ###### ######
  rstan::rstan_options(auto_write = TRUE) # save the compiled executable to '.'
  options(mc.cores = parallel::detectCores())

  ###### ###### RUN MODEL ###### ###### ######

# fit_stan <-
  rstan::stan(
    file = system.file("inst/rstan/covid_stan_script_MHC_V2.stan", mustWork=TRUE),
    control = list(adapt_delta = 0.92, max_treedepth = 12),
    data = defaultData(),
    seed = 1234,
    chains = 3,
    iter = 500,
    warmup = 400
  )
}

# samps <- rstan::extract(fit_stan)
# 
# dev.off()
# 
# system(paste("open", pdfnam))
