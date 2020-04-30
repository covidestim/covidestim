#' Configure a Covidcast run on a set of data and priors
#'
#' This function sets up the model on a set of data and priors
#'
#' @param chains The number of chains to use
#' @param iter The number of iterations to run
#' @param N_days A number. The number of days of data being modeled.
#' @param N_days_delay. A number. How many days before the first day of model
#'   data should be modeled?
#' @param seed A number. The random number generator seed for use in sampling.
#'
#' @return An S3 object of type \code{covidcast}. This can be passed to 
#'   \code{\link{run}} to execute the model. This object can also be saved
#'   to disk using \code{\link[base]{saveRDS}} to enable model replicability.
#'
#' @examples
#' covidcast(N_days = 50, seed = 42)
#' @importFrom magrittr %>%
#' @export
covidcast <- function(chains=3, iter=500, N_days, N_days_delay=10, seed=1234) {

  att(is.numeric(N_days), N_days >= 1)


  config <- defaultConfig(N_days=N_days, N_days_delay=N_days_delay)
  config <- splice_class(config, list(N_days=N_days), 'modelconfig')

  # All user-specified config-related things must be specified above this line
  # to avoid double-validation/no-validation
  if (!missing(N_days))
    validate.modelconfig(config)

  list(
    config  = config,
    chains  = chains,
    iter    = iter,
    warmup  = round(0.8*iter), # Warmup runs should be 80% of iter runs
    file    = "stan_program_default.stan",
    seed    = seed,
    control = list(adapt_delta = 0.92, max_treedepth = 12)
  ) -> properties

  structure(properties, class='covidcast')
}

#' @export
#' @rdname run.covidcast
run <- function(...) UseMethod('run')

#' @export
run.default <- function(...) stop("Must pass an object of type `covidcast`")

#' Run the Covidcast model
#'
#' Calling \code{run()} with a \code{\link{covidcast}} object executes the
#' model and returns a result. The first time the model is run, the C++
#' executable will need to be compiled, and model sampling will not begin
#' until compilation is done (typically 1-5 minutes). Afterwards, a cached
#' copy of the executable can be used, speeding execution. \code{run.covidcast}
#' will attempt to run on as many cores as appear to be available on the host
#' machine, through calling \code{\link[parallel]{detectCores}}.
#'
#' @param cc A valid \code{\link{covidcast}} configuration
#'
#' @return A S3 object of class \code{covidcast_result} containing the
#'   configuration used to run the model, the raw results, and the extracted
#'   result as produced by \code{\link[rstan]{extract}}.
#'
#' @export
run.covidcast <- function(cc) {
  # save the compiled executable to '.'
  rstan::rstan_options(auto_write = TRUE) 
  options(mc.cores = parallel::detectCores())

  rstan::stan(
    file    = system.file(paste0("rstan/", cc$file),
                          package='covidcast',
                          mustWork=TRUE),
    control = cc$control,
    data    = cc$config,
    seed    = cc$seed,
    chains  = cc$chains,
    iter    = cc$iter,
    warmup  = cc$warmup
  ) -> result

  fitted <- rstan::extract(result)

  structure(
    list(result    = result,
         summary   = rstan::summary(result)$summary,
         extracted = rstan::extract(result),
         config    = cc$config),
    class='covidcast_result'
  )
}

#' @export
"+.covidcast" <- function(a, b) {
  # 'a' is covidcast, 'b' should be priors or modelconfig
  covidcast_add(b, a)
}

covidcast_add <- function(rightside, leftside) UseMethod('covidcast_add')

#' When adding priors, we want to be sure that a new 'modelconfig' object is
#' created, in order to check these priors
#' @importFrom glue glue
covidcast_add.priors <- function(rightside, leftside) {
  newConfig       <- leftside$config + rightside
  leftside$config <- newConfig
  structure(leftside, class='covidcast')
}

covidcast_add.input <- function(rightside, leftside) {
  newConfig       <- leftside$config + rightside
  leftside$config <- newConfig
  structure(leftside, class='covidcast')
}

#' @export
print.covidcast <- function(cc) {
'Covidcast Configuration:

Seed:\t{cc$seed}
Chains:\t{cc$chains}
Iterations:\t{cc$iter}
Warmup runs:\t{cc$warmup}
Priors: Valid
Stan file:\t{cc$file}
N_days:\t{cc$config$N_days}


' -> model_summary

  substituted_string <- glue(model_summary)

  cat(substituted_string)

  # Prepend priors/inputs with a tab so that they get indented
  print.priors(cc$config, .tab = TRUE)
  print.inputs(cc$config, .tab = TRUE)
}
