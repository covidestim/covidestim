#' The 'covidcast' package.
#'
#' @description A DESCRIPTION OF THE PACKAGE
#'
#' @docType package
#' @name covidcast
#' @useDynLib covidcast, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#'
#' @references
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.19.3. https://mc-stan.org
#'
NULL

#' Configure a Covidcast run on a set of data and priors
#'
#' \code{covidcast} returns a base configuration of the model with the default
#' set of priors, and no input data. This configuration, after adding input
#' data (see \code{\link{input_cases}}, \code{\link{input_deaths}}, and
#' \code{\link{input_fracpos}}, represents a valid model configuration that can
#' be passed to \code{\link{run}}.
#'
#' @param chains The number of chains to use during MCMC, as passed to
#'   \code{\link[rstan]{sampling}}.
#' @param iter The number of iterations to run during MCMC, as passed to
#'   \code{\link[rstan]{sampling}}.
#' @param thin A positive integer to specify period for saving samples, as
#'   passed to \code{\link[rstan]{sampling}}. 
#' @param N_days A positive integer. The number of days of input data being
#'   modeled. This should always be set to the number of days in your input data.
#' @param N_days_before A positive integer. How many days before the first day
#'   of model data should be modeled?
#' @param rho_sym A number in \code{(0, 1]}. Needs documentation.
#' @param rho_sev A number in \code{(0, 1]}. Needs documentation.
#' @param seed A number. The random number generator seed for use in sampling.
#'
#' @return An S3 object of type \code{covidcast}. This can be passed to 
#'   \code{\link{run}} to execute the model. This object can also be saved
#'   to disk using \code{\link[base]{saveRDS}} to enable reproducibility across
#'  platforms or sessions.
#'
#' @examples
#' covidcast(N_days = 50, seed = 42)
#' @importFrom magrittr %>%
#' @export
covidcast <- function(N_days, N_days_before=28,
                      chains=3, iter=1500, thin = 1, 
                      rho_sym = 1, rho_sev = 0.5, seed=42) {

  att(is.numeric(N_days), N_days >= 1)

  defaultConfig(
    N_days = N_days,
    N_days_before = N_days_before,
    rho_sym = rho_sym, 
    rho_sev = rho_sev
  ) -> config

  # All user-specified config-related things must be specified above this line
  # to avoid double-validation/no-validation
  if (!missing(N_days) || !missing(N_days_before))
    validate.modelconfig(config)

  list(
    config  = config,
    chains  = chains,
    iter    = iter,
    thin    = thin,
    warmup  = round(0.8*iter), # Warmup runs should be 80% of iter runs
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
#' model and returns a result. \code{run} will attempt to run on as
#' many cores as appear to be available on the host machine, through calling
#' \code{\link[parallel]{detectCores}}. Model runtimes will range anywhere from
#' 5-120 minutes.
#'
#' When running in an interactive/TTY environment (like Rstudio, Radian, or the
#' R terminal), progress messages from \code{rstan} will be displayed,
#' indicatng how many iterations each chain has completed. When running in
#' other environments, for instance on a cluster, \code{rstan} will produce 
#' no output until the end of sampling.
#'
#' @param cc A valid \code{\link{covidcast}} configuration
#' @param cores A number. How many cores to use to execute runs.
#' @param ... Extra arguments to be passed to \code{\link[rstan]{sampling}}
#'
#' @return A S3 object of class \code{covidcast_result} containing the
#'   configuration used to run the model, the raw results, the extracted
#'   result as produced by \code{\link[rstan]{extract}}, and the summarized
#'   results as produced by \code{\link[rstan]{extract}}.
#'
#' @export
run.covidcast <- function(cc, cores = parallel::detectCores(), ...) {

  if (is.null(cc$config$obs_cas))
    stop("Case data was not entered. See `?input_cases`.")

  if (is.null(cc$config$obs_die))
    stop("Deaths data was not entered. See `?input_deaths`.")
  
  rstan::sampling(
    object  = stanmodels$stan_program_default,
    data    = structure(cc$config, class="modelconfig"),
    cores   = cores,
    control = cc$control,
    seed    = cc$seed,
    chains  = cc$chains,
    iter    = cc$iter,
    thin    = cc$thin,
    warmup  = cc$warmup,
    ...
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
  # 'a' is covidcast, 'b' should be priors or input
  covidcast_add(b, a)
}

covidcast_add <- function(rightside, leftside) UseMethod('covidcast_add')

#' When adding priors, we want to be sure that a new 'modelconfig' object is
#' created, in order to check these priors
#' @importFrom glue glue
covidcast_add.priors <- function(rightside, leftside) {
  newConfig       <- structure(leftside$config, class="modelconfig") + rightside
  leftside$config <- newConfig
  structure(leftside, class='covidcast')
}

covidcast_add.input <- function(rightside, leftside) {
  newConfig       <- structure(leftside$config, class="modelconfig") + rightside
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

  print.modelconfig(cc)
}
