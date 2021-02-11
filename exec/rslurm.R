library(rslurm)
library(covidestim)

library(docopt)
library(readr)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(assertthat)
library(purrr, warn.conflicts = FALSE)
library(glue,  warn.conflicts = FALSE)
library(cli)
library(stringr)

options(warn = 1)

glue('covidestim SLURM batch submitter

Usage:
  {name} [--cpus-per-task=<cores>] [--id-vars=<vars>] [--partition=<partition>] (--time=<time>) (--name=<jobname>) <path>
  {name} (-h | --help)
  {name} --version

Options:
  -h --help                 Show this screen.
  --cpus-per-task=<cores>   How many cores should each task (process) use?
  --id-vars=<vars>          Grouping vars in <path>, separated by commas. [default: state]
  --partition=<partition>   Which partition to run on [default: covid]
  --time=<time>             How long to run the job for
  --name=<jobname>          Name of SLURM job
  --version                 Show version.
', name = "MPIrun.R") -> doc

arguments <- docopt(doc, version = 'covidestim SLURM batch submitter 0.1')
# print(arguments)

input_file_path  <- arguments$path
jobname          <- arguments$name
partition        <- arguments$partition
id_vars          <- str_split(arguments$id_vars, ',')[[1]]
cpus_per_task    <- as.numeric(arguments$cpus_per_task)
runtime          <- arguments$time # Stored as a string
ntasks           <- as.numeric(arguments$ntasks)

# Read input data, make sure you capture the date properly, and be sure that
# there aren't any missing values in the cases/deaths data
read_csv(
  input_file_path,
  col_types = cols(date = col_date(format = "%Y-%m-%d"))
) -> d

assert_that(!any(is.na(d$cases)))
assert_that(!any(is.na(d$deaths)))

cli_alert_success("Read {.file {input_file_path}}, all assertions passed")

# This function takes some data for a geographic tract (likely a state), and
# produces a `covidestim` object which, if created without warnings or errors,
# represents a valid model configuration that is ready to be run.
#
# `type_cases`, `type_deaths` should be either "reported" or "occurred"
gen_covidestim_run <- function(d, type_cases = "reported",
                              type_deaths = "reported", ndays_before = 28) {

  d_cases   <- transmute(d, date, observation = cases)
  d_deaths  <- transmute(d, date, observation = deaths)
  ndays     <- nrow(d_cases)

  covidestim::covidestim(ndays = ndays,
                         ndays_before = ndays_before) +
    covidestim::input_cases(d_cases,   type = type_cases) +
    covidestim::input_deaths(d_deaths, type = type_deaths) +
}

# Take the data, and split it into smaller tibbles for each state
runs_nested <- group_by_at(d, id_vars) %>% nest() %>%
  mutate(config = map(data, gen_covidestim_run),
         group_name = paste(as.list(cur_group()), collapse = '-'))

cli_alert_success("Configuration generated successfully. Summary:")
cli_ul()
pwalk(runs_nested, function(data, config, group_name, ...) {
  cli_li(group_name)
  ulid <- cli_ul()
  cli_li("Death data: {nrow(data)} obs")
  cli_li("Cases data: {nrow(data)} obs")
  cli_end(ulid)
})
cli_end()

do_covidestim_run <- function(config, group_name) {
  safe_run           <- purrr::quietly(covidestim::run) # capture warning msgs
  diagnostic_file    <- glue::glue("logs/stan-{group_name}") # Stan logs
  startTime          <- Sys.time()

  cli::cli_alert_info("Started {group_name} at {startTime}")

  safe_run(
    config,
    cores = cpus_per_task,
    diagnostic_file = diagnostic_file
  ) -> result

  endTime <- Sys.time()
  cli::cli_alert_info("Finished {group_name} on {core} at {endTime}")
  cli::cli_alert_info("Runtime was {prettyunits::pretty_dt(endTime - startTime)}")

  warnings_observed <- result$warnings

  tibble::tibble(group_name,
                 result   = list(result),
                 time     = endTime - startTime,
                 warnings = warnings_observed)
}

print(jobname)

slurm_apply(
  do_covidestim_run,
  runs_nested %>% ungroup() %>% select(config, group_name),
  add_objects   = "cpus_per_task",
  jobname       = jobname,
  nodes         = nrow(runs_nested),
  cpus_per_node = 1, # HACK: Force rslurm to not do weird multithreading stuff,
                     # But then repeat the argument later with the real 
                     # `cpus_per_task` value, which is the one that SLURM will
                     # actually look at, since rightmost args have precedence
  preschedule_cores = FALSE,
  sh_template   = "batch_template.txt",
  slurm_options = list(
    `cpus-per-task` = cpus_per_task,
    `mem-per-cpu`   = "1G",
    time            = runtime,
    partition       = partition
  ),
  submit = FALSE
)

cli_alert_success("Finished model setup\n")
