#' Run Covidcast on a set of data and priors
#'
#' This function executes the model on a set of data and priors
#'
#' @param config Description
#'
#' @param chains Description
#'
#' @param iter Description
#'
#' @return The return value
#'
#' @examples
#' print(mtcars)
#' @importFrom magrittr %>%
#' @export
covidcast <- function(config=defaultConfig(), chains=3, iter=500) {

  # save the compiled executable to '.'
  rstan::rstan_options(auto_write = TRUE) 
  options(mc.cores = parallel::detectCores())

  result <- rstan::stan(
    file    = system.file("rstan/covid_stan_script_MHC_V2.stan",
                          package="covidcast",
                          mustWork=TRUE),
    control = list(adapt_delta = 0.92, max_treedepth = 12),
    data    = config,
    seed    = 1234,
    chains  = chains,
    iter    = iter,
    warmup  = round(0.8*iter) # Warmup runs should be 80% of iter runs
  )

  fitted <- rstan::extract(result)

  structure(
    list(result=result,
         fitted=fitted,
         config=config),
    class='covidcast'
  )
}
