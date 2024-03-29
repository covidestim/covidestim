#' The 'covidestim' package.
#'
#' @description Real-time estimates of the true size and trajectory of local
#'   COVID-19 epidemics.
#'
#' @docType package
#' @name covidestim-package
#' @useDynLib covidestim, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom Rdpack reprompt
#'
#' @references
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.19.3. https://mc-stan.org
#'
NULL

#' Configure a Covidestim run on a set of priors and input data
#'
#' \code{covidestim} returns a base configuration of the model with the default
#' set of priors, and no input data. This configuration, after adding input
#' data (see \code{\link{input_cases}}, \code{\link{input_deaths}},
#' \code{priors_*}), represents a valid model configuration that can be passed
#' to \code{\link{run.covidestim}} (for NUTS) or
#' \code{\link{runOptimizer.covidestim}} (for BFGS).
#'
#' @param chains The number of chains to use during NUTS sampling, as passed to
#'   \code{\link[rstan]{sampling}}.
#' @param iter The number of iterations to run during NUTS sampling, as passed
#'   to \code{\link[rstan]{sampling}}.
#' @param thin A positive integer to specify period for saving samples, as
#'   passed to \code{\link[rstan]{sampling}}. Modify this only if you intend
#'   to inspect raw iterations data returned by \code{\link[rstan]{sampling}}.
#' @param nweeks A positive integer. The number of weeks of input data being
#'   modeled. This should always be set to the number of weeks in your input
#'   data.
#' @param nweeks_before A positive integer. How many weeks before the first week
#'   of model data should be modeled? A higher number will produce estimates
#'   that go farther back in time, however, those estimates will contain more
#'   and more uncertainty.
#' @param pop_size A positive integer What is the population in the geography 
#'   being modelled? This sets the max susceptible population and becomes
#'   important as the population ever infected approaches the population size.
#' @param seed A number. The random number generator seed for use in NUTS
#'   sampling or in BFGS.
#' @param region A string. The FIPS code (for U.S. counties) or state name
#'   (e.g. \code{New York}) being modeled. Required.
#' @param nspl_rt_knotwidth An integer. The knotwidth to use for the spline of 
#' Rt.
#' @param start_p_imm A probability (number between 0 and 1). The starting 
#'    fraction of the population with immunity against the Omicron variant.
#'    Required.
#' @param cum_p_inf_init A probability (number between o and 1). The starting
#'    fraction of the population that has ever been infected. Required.
#'
#' @return An S3 object of type \code{covidestim}. This can be passed to
#' \code{\link{run.covidestim}} or \code{\link{runOptimizer.covidestim}} to execute the model, as
#' long as input data has been added (using the addition operator, see
#' example). This object can also be saved to disk using
#' \code{\link[base]{saveRDS}} to enable reproducibility across platforms or
#' sessions. However, note that Stan runs are only reproducible under very
#' specific conditions due to Stan's multi-language architecture. The
#' \code{print} method is overloaded to return to the user a summary of the
#' configuration, including prior values and the presence or absence of input
#' data.
#'
#' @seealso \url{https://mc-stan.org/docs/2_18/reference-manual/reproducibility-chapter.html}
#'
#' @examples
#' cfg <- covidestim(nweeks = 31, region = 'Connecticut',
#'    pop = get_pop("Connecticut"),
#'    start_p_imm = get_imm_init("Connecticut")$start_p_imm,
#'    cum_p_inf_init = get_imm_init("Connecticut")$cum_p_inf_init) +
#'   input_cases(example_ct_data('cases')) +
#'   input_deaths(example_ct_data('deaths')) +
#'   input_rr(example_ct_data('RR')) + 
#'   input_hosp(example_ct_data('hosp')) +
#'   input_boost(example_ct_data('boost'))
#'   
#'
#' @importFrom magrittr %>%
#' @export
covidestim <- function(nweeks, 
                       nweeks_before = 4,
                       pop_size = 1e12,
                       chains = 3, 
                       iter = 2000, 
                       nspl_rt_knotwidth = 1,
                       thin = 1, 
                       seed = 42,
                       adapt_delta = 0.98, 
                       max_treedepth = 14,
                       window.length = 7,
                       region,
                       start_p_imm,
                       cum_p_inf_init) {

  att(is.numeric(nweeks), nweeks >= 1)

  defaultConfig(
    N_weeks            = nweeks,
    N_weeks_before     = nweeks_before,
    pop_size           = pop_size,
    n_spl_rt_knotwidth = nspl_rt_knotwidth,
    region             = region,
    start_p_imm        = start_p_imm,
    cum_p_inf_init     = cum_p_inf_init
  ) -> config

  # All user-specified config-related things must be specified above this line
  # to avoid double-validation/no-validation
  # if (!missing(ndays) || !missing(ndays_before))
  #   validate.modelconfig(config)
  if (!missing(nweeks) || !missing(nweeks_before))
    validate.modelconfig(config)

  list(
    config  = config,
    chains  = chains,
    iter    = iter,
    thin    = thin,
    warmup  = round((2/3)*iter), # Warmup runs should be 66% of iter runs
    seed    = seed,
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)
  ) -> properties

  structure(properties, class='covidestim')
}

#' Population estimates for US states and counties
#'
#' Returns 2019 census estimate of state or county population.
#'
#' @param region A string with the state name, or the FIPS code
#'
#' @return State/county population as a numeric, or an error
#'
#' @examples
#' get_pop('Connecticut')
#' get_pop('09009')
#'
#' @export
get_pop <- function(region) {
  found <- dplyr::filter(pop_state, state == region)

  if (nrow(found) == 0)
    found <- dplyr::filter(pop_county, fips == region)

  if (nrow(found) == 0)
    stop(glue::glue("Could not find population data for region {region}!"))

  if (nrow(found) > 1)
    stop(glue::glue("Found more than set of population data for region {region}!"))

  found$pop
}

#' Immunity and cumulative infections (fraction) estimates for US states and counties 
#'
#' Returns December 1, 2021 covidestim estimates of state or county effective protection
#' and fraction of the population ever infected.
#'
#' @param region A string with the state name, or the FIPS code
#'
#' @return State/county population as a numeric, or an error
#'
#' @examples
#' get_imm_init('Connecticut')
#' get_imm_init('09009')
#'
#' @export
get_imm_init <- function(region) {
  found <- dplyr::filter(imm_state, location == region)

  if (nrow(found) == 0)
    found <- dplyr::filter(imm_county, location == region)

  if (nrow(found) == 0)
    stop(glue::glue("Could not find immunity data for region {region}!"))

  if (nrow(found) > 1)
    stop(glue::glue("Found more than set of immunity data for region {region}!"))

  found
}


#' @export
run <- function(...) UseMethod('run')

#' Run the Covidestim model using NUTS
#'
#' Calling \code{run()} with a \code{\link{covidestim}} object executes the
#' model and returns a result. \code{run} will attempt to run on as
#' many cores as appear to be available on the host machine, through calling
#' \code{\link[parallel]{detectCores}}. Model runtimes will range anywhere from
#' 30 minutes to 12 hours, with the cumulative number of cases appearing to be
#' a strong correlate to longer runtime.
#'
#' When running in an interactive/TTY environment (like Rstudio, Radian, or the
#' R terminal), progress messages from \code{rstan} will be displayed,
#' indicating how many iterations each chain has completed. When running in
#' other environments, for instance on a cluster, \code{rstan} will produce 
#' no output until the end of sampling.
#'
#' The sampler may return warnings upon completion. In general, "treedepth" and
#' "divergent transitions" are the most serious. Use caution when interpreting
#' results which were accompanied by these messages. See [mc-stan.org](https://mc-stan.org/misc/warnings.html) for
#' a detailed guide to Stan's most common warning messages.
#'
#' A second method for fitting the model, using the BFGS algorithm, is available
#' as [runOptimizer]. It is significantly faster, but lacks CI's.
#'
#' | **Function** | **Method** | **CI's** | **Speed** | **Termination** |
#' | ---      | ---         | ---  | ---   | ---         |
#' | `run()`  | NUTS        | Yes  | 30m-hours | Always, potentially with warnings, of which "treedepth" and "divergent transitions" are the most serious |
#' | `runOptimizer()` | BFGS | No | ~1-3min | Potentially with nonzero exit status (lack of convergence), or timeout (rare, gracefully handled internally)
#'
#' @param cc A valid \code{\link{covidestim}} configuration
#' @param cores A number. How many cores to use to execute runs.
#' @param ... Extra arguments, passed on to \code{\link[rstan]{sampling}}
#'
#' @return An S3 object of class \code{covidestim_result} containing the
#'   configuration used to run the model, the raw Stan results, the extracted
#'   result as produced by \code{\link[rstan]{extract}}, and the summarized
#'   results as produced by \code{\link[rstan]{extract}}.
#'
#' @seealso [runOptimizer.covidestim]
#'
#' @examples
#' cfg <- covidestim(nweeks = 31, region = 'Connecticut',
#'    pop = get_pop("Connecticut"),
#'    start_p_imm = get_imm_init("Connecticut")$start_p_imm,
#'    cum_p_inf_init = get_imm_init("Connecticut")$cum_p_inf_init) +
#'   input_cases(example_ct_data('cases')) +
#'   input_deaths(example_ct_data('deaths')) +
#'   input_rr(example_ct_data('RR')) + 
#'   input_hosp(example_ct_data('hosp')) +
#'   input_boost(example_ct_data('boost'))
#'   
#' 
#' \dontrun{
#' result <- run(cfg, cores = 2)
#' }
#'
#' @export
run.covidestim <- function(cc, cores = parallel::detectCores(), ...) {

    # Require that case and death data be entered
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
         config    = cc$config,
         flags     = ""),
    class='covidestim_result'
  )
}

#' @export
runOptimizer <- function(...) UseMethod('runOptimizer')

#' Run the Covidestim model using BFGS
#'
#' In addition to NUTS sampling, you can fit the model using the BFGS algorithm.
#' This algorithm tries to maximize the mode of the posterior; it runs much
#' faster than NUTS, but won't return any confidence intervals. On
#' \url{https://covidestim.org}, BFGS is used to produce estimates for counties, and as
#' a fallback when NUTS fails to produce a timely and/or converged fit for state
#' data. Runtimes scale linearly to the value of \code{tries/cores}, but are
#' generally in the 30s-10min range. The same underlying model is used.
#'
#' The BFGS algorithm is run \code{tries} times using \code{tries} different
#' seeds derived from \code{cc$seed}. Once all runs complete, the run with the
#' highest value of the log-posterior is selected and returned. BFGS
#' occasionally gets stuck; these runs are flagged and discarded before the
#' ranking begins.
#'
#' A second method for fitting the model, using the NUTS algorithm, is available
#' as [run]. It provides CI's, but is significantly slower.
#'
#' | **Function** | **Method** | **CI's** | **Speed** | **Termination** |
#' | ---      | ---         | ---  | ---   | ---         |
#' | `run()`  | NUTS        | Yes  | 30m-hours | Always, potentially with warnings, of which "treedepth" and "divergent transitions" are the most serious |
#' | `runOptimizer()` | BFGS | No | ~1-3min | Potentially with nonzero exit status (lack of convergence), or timeout (rare, gracefully handled internally)
#'
#' @param cc A valid \code{\link{covidestim}} configuration
#' @param cores A number. How many cores to use to execute runs. Multicore
#'   execution is less important for \code{runOptimizer} than it is for
#'   \code{\link{run}}.
#' @param tries The number of times to run the BFGS algorithm.
#' @param iter Passed to \code{\link[rstan]{optimizing}}.
#' @param timeout How long to let each run go for before killing it, in seconds.
#' @param ... Arguments forwarded to \code{\link[rstan]{optimizing}}.
#'
#' @return An S3 object of class \code{covidestim_result}
#'
#' @seealso [run.covidestim]
#'
#' @examples
#' cfg <- covidestim(nweeks = 31, region = 'Connecticut',
#'    pop = get_pop("Connecticut"),
#'    start_p_imm = get_imm_init("Connecticut")$start_p_imm,
#'    cum_p_inf_init = get_imm_init("Connecticut")$cum_p_inf_init) +
#'   input_cases(example_ct_data('cases')) +
#'   input_deaths(example_ct_data('deaths')) +
#'   input_rr(example_ct_data('RR')) + 
#'   input_hosp(example_ct_data('hosp')) +
#'   input_boost(example_ct_data('boost'))
#' 
#' \dontrun{
#' result <- runOptimizer(cfg, cores = 2)
#' }
#'
#' @export
runOptimizer.covidestim <- function(cc,
                                    cores   = parallel::detectCores(),
                                    tries   = 10,
                                    iter    = 2000,
                                    timeout = 60,
                                    ...) {
  
  # Require that case and death data be entered
  if (is.null(cc$config$obs_cas))
    stop("Case data was not entered. See `?input_cases`.")
  
  if (is.null(cc$config$obs_die))
    stop("Deaths data was not entered. See `?input_deaths`.")
  
  # Set the RNG seed to the seed specified in the `covidestim_config` object.
  # Then, use that seed to create `tries` random integers, each of which will
  # be used as the seed value for an instance of `rstan::optimizing`.
  set.seed(cc$seed)
  seeds <- sample.int(.Machine$integer.max, tries)

  runOptimizerWithSeed <- function(seed, i) {
    startTime <- Sys.time()

    rstan::optimizing(
      object    = stanmodels$stan_program_default,
      data      = structure(cc$config, class="modelconfig"),
      algorithm = "BFGS",
      seed      = seed,
      iter      = iter,
      tol_obj   = 1e-10,
      as_vector = FALSE, # Otherwise you get a sloppy list structure
      ...
    ) -> result

    endTime <- Sys.time()

    message(glue::glue(
      'Finished try #{i} in {dt} with exit code {ec}',
      dt = prettyunits::pretty_dt(endTime - startTime),
      ec = result$return_code
    ));

    result
  }

  # This function will return NULL when there is a timeout
  runOptimizerWithSeedInTime <- function(timeout, ...)
    tryCatch(
      R.utils::withTimeout(runOptimizerWithSeed(...), timeout = timeout),
      error = function(c) {
        message(glue::glue(
          'Abandoned try #{i} due to timeout',
          i = list(...)[[2]]
        ))

        NULL
      }
    )

  # By default, use a sequential mapping function
  map_fun <- purrr::imap

  # Set up multicore execution for the `furrr::future_imap` call.
  # Replace the mapping function with a parallel mapping function that has the
  # same functional form, but executes in parallel.
  if (cores > 1) {
    future::plan(future::multisession, workers = cores)
    map_fun <- furrr::future_imap
  }

  # Run the optimizer on all the different seeds
  results <- map_fun(
    seeds,
    runOptimizerWithSeedInTime,
    timeout = timeout
  )

  # Change the execution plan back to sequential
  if (cores > 1)
    future::plan(future::sequential)

  # Return code of 0 indicates success for `rstan::optimizing`. This is just
  # a standard UNIX return code b/c `rstan::optimizing` calls into CmdStan.
  successful_results <-
    purrr::discard(results, is.null) %>% # Removes timed-out runs
    purrr::keep(., ~.$return_code == 0)  # Removes >0 return-val runs

  if (length(successful_results) == 0) {
    rlang::abort(
      "all_runs_were_bad_error",
      message = "All BFGS runs timed out or failed!",
      path = "R/covidestim.R"
    )
  }

  # Extract the mode of the posterior from the results that didn't time out
  # and didn't return an error code of 70
  opt_vals <- purrr::map_dbl(successful_results, 'value') 

  # In theory the log posterior could be infinite (likely, -Infinity), which
  # wouldn't be valid but would technically be the maximum value. Throw an
  # error in this case.
  if (is.infinite(max(opt_vals)))
    stop(glue::glue(
      'The value of the log posterior was infinite for these runs:\n{runs}',
      runs = which(is.infinite(opt_vals) & opt_vals > 0)
    ))

  # The first successful result which has `opt_val` equal to the maximum
  # `opt_val` is the result that will be returned too the user. Note that it's
  # unlikely that there will be more that one trajectory with the same
  # `opt_val`. However, if this is the case, the first of these results will
  # be returned
  result <- successful_results[which(opt_vals == max(opt_vals))][[1]]

  c(
    "deaths",
    "deaths_of_diagnosed",
    "diagnoses",
    "diagnoses_of_symptomatic",
    "effective_protection_inf_prvl",
    "effective_protection_inf_vax_boost_prvl",
    "effective_protection_inf_vax_prvl",
    "effective_protection_vax_boost_prvl",
    "effective_protection_vax_prvl",
    "fitted_cases",
    "fitted_deaths",
    "fitted_hospitalizations",
    "fitted_wastewater_prvl",
    "immunoexposed_cumulative",
    "infections_cumulative",
    "infections",
    "infections_premiere",
    "r_t",
    "seropositive_prvl",
    "severe",
    "susceptible_prvl",
    "susceptible_severe_prvl",
    "symptomatic"
  ) -> essential_vars

  # c(
  #   "p_diag_if_sym",
  #   "p_diag_if_asy",
  #   "p_diag_if_sev",
  #   "p_sym_if_inf",
  #   "p_sev_if_sym"
  # ) -> extraParams


  structure(
    #list(result   = result$par[c(essential_vars, extraParams)],
    list(result   = result$par[essential_vars],
         opt_vals = opt_vals,
         config   = cc$config,
         flags    = "optimizer"),
    class='covidestim_result'
  )
}

#' @export
"+.covidestim" <- function(a, b) {
  # 'a' is covidestim, 'b' should be priors or input
  covidestim_add(b, a)
}

covidestim_add <- function(rightside, leftside) UseMethod('covidestim_add')

#' When adding priors, we want to be sure that a new 'modelconfig' object is
#' created, in order to check these priors
covidestim_add.priors <- function(rightside, leftside) {
  newConfig       <- structure(leftside$config, class="modelconfig") + rightside
  leftside$config <- newConfig
  structure(leftside, class='covidestim')
}

covidestim_add.input <- function(rightside, leftside) {
  newConfig       <- structure(leftside$config, class="modelconfig") + rightside
  leftside$config <- newConfig
  structure(leftside, class='covidestim')
}

#' @export
print.covidestim <- function(cc) {
'Covidestim Configuration:

Seed:\t{cc$seed}
Chains:\t{cc$chains}
Iterations:\t{cc$iter}
Warmup runs:\t{cc$warmup}
Priors: Valid
nweeks:\t{cc$config$N_weeks}


' -> model_summary

  substituted_string <- glue(model_summary)

  cat(substituted_string)

  print.modelconfig(cc)
}
