###### ###### SETTINGS ###### ###### ######
rstan_options(auto_write = TRUE) # save the compiled executable to '.'
options(mc.cores = parallel::detectCores())

###### ###### RUN MODEL ###### ###### ######

fit_stan <-
  stan(file = "covid_stan_script_MHC_V2.stan",
       control = list(adapt_delta = 0.92, max_treedepth = 12),
       data = defaultData(),
       seed = 1234,
       chains = 3,
       iter = 500,
       warmup = 400)

samps <- extract(fit_stan)

# summary(fit_stan) names(samps)


dev.off()

system(paste("open", pdfnam))

######################## 
