#!/usr/bin/env Rscript
'covidestim-batch

Usage:
  covidestim-batch --input <path> --key <key> --metadata <path> [--attempts <num>] [--always-sample] [--always-optimize] [--cpus <num>] --save-summary <path> [--save-warning <path>] [--save-optvals <path>] [--save-method <path>] [--save-metadata <path>] [--save-raw <path>] [--dbstan-insert]
  covidestim-batch (-h | --help)
  covidestim-batch --version

Options:
  --input <path>          Path to CSV of input data
  --key <key>             "state" or "fips"
  --metadata <path>       Path to JSON of metadata
  --attempts <num>        Number of attempts, total [default: 1]
  --always-sample         Always use the sampler.
  --always-optimize       Always use the optimizer. Implies attempt=1
  --cpus <num>            How many CPUs to use [default: 1]
  --save-summary <path>   Where to save a CSV of the summary dataframe
  --save-warning <path>   Where to save a CSV of all warnings
  --save-optvals <path>   Where to save a CSV of all optvals
  --save-method <path>    Where to save a CSV of which method was ultimately used
  --save-metadata <path>  Where to save a JSON of metadata
  --save-raw <path>       Where to save an RDS archive of each result
  --dbstan-insert         Insert sampled `stanfit` objects into a dbstan database
  -h --help  Show this screen.
  --version  Show version.

The following environment variables must be set to use `--dbstan-insert`:
- COVIDESTIM_DBSTAN_HOST
- COVIDESTIM_DBSTAN_USER
- COVIDESTIM_DBSTAN_PASS
- COVIDESTIM_DBSTAN_DBNAME
' -> doc

suppressPackageStartupMessages({
  library(docopt)
  library(tidyverse)
  library(jsonlite)
  library(cli)
  library(dbstan)
  library(callr)
})
library(covidestim)

args <- docopt(doc, version = 'covidestim-batch 0.1')

runner_sampler   <- quietly(covidestim::run)
runner_optimizer <- quietly(covidestim::runOptimizer)

output_alert <- function(name, fname)
  cli_alert_info("Writing {name} to {.file {fname}}")

append_in <- function(lst, where, key, val)
  purrr::modify_in(lst, where, function(x) { x[[key]]<-val; x})

conn <- NULL
if (args$dbstan_insert)
  conn <- DBI::dbConnect(
    RPostgres::Postgres(),
    host     = Sys.getenv("COVIDESTIM_DBSTAN_HOST"),
    user     = Sys.getenv("COVIDESTIM_DBSTAN_USER"),
    password = Sys.getenv("COVIDESTIM_DBSTAN_PASS"),
    dbname   = Sys.getenv("COVIDESTIM_DBSTAN_DBNAME")
  )

input_df_colspec <- rlang::list2(
  !!args$key := col_character(),
  date        = col_date(),
  cases       = col_number(),
  deaths      = col_number(),
  RR          = col_number(),
  hosp        = col_number(),
  boost       = col_number(),
  missing_hosp= col_logical()
) %>% do.call(cols, .)

cli_alert_info("Reading input data from {.file {args$input}}")
d <- read_csv(args$input, col_types = input_df_colspec) %>%
  group_by(across(all_of(args$key)))

cli_alert_info("Reading metadata from {.file {args$metadata}}")
# Read metadata into a list, keyed on `args$key` ("state" or "fips")
metadata <- read_json(args$metadata) %>% setNames(., map_chr(., args$key))

cli_alert_info("Tracts in this process:")
cli_ul(d[, args$key] %>% unique)

aggregated_results <- group_map(d, function(input_data, group_keys) {

  region <- group_keys[[args$key]]
  cli_h1("{region}")

  get_input <- function(col_name) select(input_data, date, observation := !!col_name)

  imminits <- get_imm_init(region)

  covidestim_config_options <- list(
    nweeks         = nrow(input_data),
    seed           = sample.int(.Machine$integer.max, 1),
    region         = region,
    cum_p_inf_init = imminits$cum_p_inf_init,
    start_p_imm    = imminits$start_p_imm,
    pop_size       = get_pop(region),
    nweeks_before  = 4
  )

  lastHospDate <- as.Date(metadata[[region]]$lastHospDate)

  cfg <- do.call(covidestim, covidestim_config_options) +
    input_cases(get_input("cases")) +
    input_deaths(get_input("deaths")) +
    input_rr(get_input("RR")) + 
    input_boost(get_input("boost")) +
    input_hosp(get_input("hosp"), lastHospDate = lastHospDate)

  tries   <- 25
  timeout <- 60
  result_optimizer <- tryCatch(
    callr::r(
      runner_optimizer,
      args = list(cfg, cores = 1, tries = tries),
      timeout = timeout,
      package = TRUE
    ),
    callr_timeout_error = function(cnd) {
      message(glue::glue("{region} has timed out after {tries} tries and {timeout} seconds. Ignoring."))
      message(cnd)
      return(NULL)
    },
    callr_error = function(cnd) {
      message(
        glue::glue("Optimizer function failed. Likely all BFGS runs of {region} failed to meet runOptimizer's successful-result conditions. Ignoring.")
      )
      message(cnd)
      return(NULL)
    }
  )

  if (is.null(result_optimizer)) {
    return(list(
      run_summary = bind_cols(!!args$key := region, tibble()),
      warnings    = bind_cols(!!args$key := region, warnings = tibble()),
      opt_vals    = bind_cols(!!args$key := region, optvals  = tibble()),
      method      = bind_cols(!!args$key := region, method   = tibble()),
      raw         = NULL
    ))
  }

  cli_alert_success("Optimizer complete. Log-posterior values:")

  run_summary <- summary(result_optimizer$result)
  warnings_optimizer <- result_optimizer$warnings
  opt_vals <- result_optimizer$result$opt_vals

  cli_ol(sort(opt_vals, decreasing = TRUE))
  
if (any((run_summary$infections / 7) > get_pop(region))) {
    
    cli_alert_warning("Optimizer did not converge on valid posterior for {region}; Discarding run from results.")
    
    return(list(
      run_summary = bind_cols(!!args$key := region, tibble()),
      warnings    = bind_cols(!!args$key := region, warnings = warnings_optimizer),
      opt_vals    = bind_cols(!!args$key := region, optvals  = opt_vals),
      method      = bind_cols(!!args$key := region, method   = "optimizer"),
      raw         = NULL
    ))
    
}

  if (
    # If it's the last attempt and we aren't always sampling, OR
    (identical(as.numeric(args$attempts), 1) && !identical(args$always_sample, TRUE)) ||

    # the --always-optimize flag was passed
    identical(args$always_optimize, TRUE)
  ) {

    cli_alert_info("Returning best optimizer result; sampler will not run")

    run_metadata <- list(optimizing = list(
      covidestim_config_options = covidestim_config_options,
      warnings = I(warnings_optimizer),
      optvals  = I(opt_vals),
      tries    = tries
    ))

    metadata <<- append_in(metadata, region, "run", run_metadata)

    return(list(
      run_summary = bind_cols(!!args$key := region, run_summary),
      warnings    = bind_cols(!!args$key := region, warnings = warnings_optimizer),
      opt_vals    = bind_cols(!!args$key := region, optvals  = opt_vals),
      method      = bind_cols(!!args$key := region, method   = "optimizer"),
      raw         = result_optimizer
    ))
  }

  start <- lubridate::now()
  cli_alert_info("Sampler starting at {start}")
  result <- runner_sampler(cfg, cores = as.numeric(args$cpus))
  end <- lubridate::now()
  cli_alert_success("Sampler ends at {end}")

  run_summary <- summary(result$result)
  warnings_sampler <- result$warnings

  dbstan_id <- NULL
  if (args$dbstan_insert) {
    dbstan_id <- stanfit_insert(result$result$result, conn, include_samples = T)
    cli_alert_success("{.code dbstan} id is {.code {dbstan_id}}")
  }

  run_metadata <- list(sampling = list(
    covidestim_config_options = covidestim_config_options,
    warnings = I(warnings_sampler),

    # Unit of `duration` is seconds
    timing = list(start = start, end = end, duration = as.numeric(end - start)),
    dbstan_id = dbstan_id
  ))

  metadata <<- append_in(metadata, region, "run", run_metadata)

  # Error on treedepth warning, or any divergent transitions warning
  # indicating >= 10 divergent transitions
  if (any(str_detect(warnings_sampler, 'treedepth')) ||
      any(str_detect(warnings_sampler, ' [0-9]{2,} divergent'))) {
    cli_alert_danger("Failing due to sampler warnings:")
    print(warnings_sampler)
    quit(status=1)
  }

  return(list(
    run_summary = bind_cols(!!args$key := region, run_summary),
    warnings    = bind_cols(!!args$key := region, warnings = warnings_sampler),
    opt_vals    =    tibble(!!args$key := region, optvals = numeric()),
    method      = bind_cols(!!args$key := region, method = "sampler"),
    raw         = result
  ))
})

cli_h1("All runs complete")

output_alert("summary", args$save_summary)
write_csv(map(aggregated_results, 'run_summary') %>% bind_rows, args$save_summary)

if (!is.null(args$save_warning)) {
  output_alert("warnings", args$save_warning)
  write_csv(map(aggregated_results, 'warnings') %>% bind_rows, args$save_warning)
}

if (!is.null(args$save_optvals)) {
  output_alert("optvals", args$save_optvals)
  write_csv(map(aggregated_results, 'opt_vals') %>% bind_rows, args$save_optvals)
}

if (!is.null(args$save_method)) {
  output_alert("method used", args$save_method)
  write_csv(map(aggregated_results, 'method') %>% bind_rows, args$save_method)
}

if (!is.null(args$save_metadata)) {
  output_alert("metadata", args$save_metadata)
  write_json(
    unname(metadata),
    args$save_metadata,
    null = 'null', auto_unbox = T
  )
}

if (!is.null(args$raw)) {
  output_alert("raw model output", args$raw)
  saveRDS(map(aggregated_results, 'raw'), args$raw)
}
