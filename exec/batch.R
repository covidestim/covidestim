#!/usr/bin/env Rscript
'covidestim-batch

Usage:
  covidestim-batch --input <path> --key <key> --metadata <path> [--attempts <num>] [--always-sample] [--always-optimize] [--cpus <num>] --save-summary <path> [--save-warning <path>] [--save-optvals <path>] [--save-method <path>] [--save-metadata <path>] [--save-raw <path>] 
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
  -h --help  Show this screen.
  --version  Show version.

' -> doc

suppressPackageStartupMessages({
  library(docopt)
  library(tidyverse)
  library(jsonlite)
  library(cli)
})
library(covidestim)

args <- docopt(doc, version = 'covidestim-batch 0.1')

runner_sampler   <- quietly(covidestim::run)
runner_optimizer <- quietly(covidestim::runOptimizer)

output_alert <- function(name, fname)
  cli_alert_info("Writing {name} to {.file {fname}}")

input_df_colspec <- rlang::list2(
  !!args$key := col_character(),
  date        = col_date(),
  .default    = col_number()
) %>% do.call(cols, .)

cli_alert_info("Reading input data from {.file {args$input}}")
d <- read_csv(args$input, col_types = input_df_colspec) %>%
  group_by(across(all_of(args$key)))

cli_alert_info("Reading metadata from {.file {args$metadata}}")
metadata <- read_json(args$metadata, simplifyVector = TRUE)

cli_alert_info("Tracts in this process:")
cli_ul(d[, args$key] %>% unique)

aggregated_results <- group_map(d, function(input_data, group_keys) {

  region <- group_keys[[args$key]]
  cli_h1("{region}")

  get_input <- function(col_name) select(input_data, date, observation := !!col_name)

  imminits <- get_imm_init(region)

  cfg <- covidestim(
    nweeks         = nrow(input_data),
    seed           = sample.int(.Machine$integer.max, 1),
    region         = region,
    cum_p_inf_init = imminits$cum_p_inf_init,
    start_p_imm    = imminits$start_p_imm,
    pop_size       = get_pop(region),
    nweeks_before  = 4
  ) +
    input_cases(get_input("cases")) +
    input_deaths(get_input("deaths")) +
    input_rr(get_input("RR")) + 
    input_boost(get_input("boost")) +
    input_hosp(get_input("hospi"))

  result_optimizer <- runner_optimizer(cfg, cores = 1, tries = 25)
  cli_alert_success("Optimizer complete. Log-posterior values:")

  run_summary <- summary(result_optimizer$result)
  warnings_optimizer <- result_optimizer$warnings
  opt_vals <- result_optimizer$result$opt_vals

  cli_ol(sort(opt_vals, decreasing = TRUE))

  if (
    # If it's the last attempt and we aren't always sampling, OR
    (identical(as.numeric(args$attempts), 1) && !identical(args$always_sample, TRUE)) ||

    # the --always-optimize flag was passed
    identical(args$always_optimize, TRUE)
  ) {

    cli_alert_info("Returning best optimizer result; sampler will not run")

    return(list(
      run_summary = bind_cols(!!args$key := region, run_summary),
      warnings    = bind_cols(!!args$key := region, warnings = warnings_optimizer),
      opt_vals    = bind_cols(!!args$key := region, optvals  = opt_vals),
      method      = bind_cols(!!args$key := region, method   = "optimizer"),
      raw         = result_optimizer
    ))
  }

  cli_alert_info("Sampler starting at {lubridate::now()}")
  result <- runner_sampler(cfg, cores = as.numeric(args$cpus))
  cli_alert_success("Sampler ends at {lubridate::now()}")

  run_summary <- summary(result$result)
  warnings_sampler <- result$warnings

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
    opt_vals    = tibble(!!args$key := region, optvals = numeric()),
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
  write_json(metadata, args$save_metadata, null = "null")
}

if (!is.null(args$raw)) {
  output_alert("raw model output", args$raw)
  saveRDS(map(aggregated_results, 'raw'), args$raw)
}
