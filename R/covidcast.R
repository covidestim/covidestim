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
covidcast <- function(config = defaultConfig(), chains=3, iter=500) {

  list(
    config  = config,
    chains  = chains,
    iter    = iter,
    warmup  = round(0.8*iter), # Warmup runs should be 80% of iter runs
    file    = "stan_program_default.stan",
    seed    = 1234,
    control = list(adapt_delta = 0.92, max_treedepth = 12)
  ) -> properties

  structure(properties, class='covidcast')
}

#' @export
run <- function(...) UseMethod('run')

#' @export
run.default <- function(...) stop("Must pass an object of type `covidcast`")

#' @export
run.covidcast <- function(cc) {
  # save the compiled executable to '.'
  rstan::rstan_options(auto_write = TRUE) 
  options(mc.cores = parallel::detectCores())

  rstan::stan(
    file    = system.file(paste0("rstan/", cc$file),
                          package=pkgload::pkg_name(),
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
         extracted = rstan::extract(result),
         config    = config),
    class='covidcast'
  )
}

#' @export
"+.covidcast" <- function(a, b) {
  # 'a' is covidcast, 'b' should be priors or modelconfig
  covidcast_add(b, a)
}

#' @export
covidcast_add <- function(rightside, leftside) UseMethod('covidcast_add')

#' When adding priors, we want to be sure that a new 'modelconfig' object is
#' created, in order to check these priors
#' @export
#' @importFrom glue glue
covidcast_add.priors <- function(rightside, leftside) {
  newConfig       <- leftside$config + rightside
  leftside$config <- newConfig
  leftside
}

print.covidcast <- function(cc) {
'Covidcast configuration

Seed:\t{cc$seed}
Chains:\t{cc$chains}
Iterations:\t{cc$iter}
Priors: Valid
Stan file:\t{cc$file}
' -> model_summary

  substituted_string <- glue(model_summary)

  cat(substituted_string)
}
