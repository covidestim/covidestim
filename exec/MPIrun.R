suppressPackageStartupMessages(library(doMPI))
library(Rmpi)
library(docopt)
library(readr)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(assertthat)
library(purrr, warn.conflicts = FALSE)
library(glue,  warn.conflicts = FALSE)
library(covidcast)
library(cli)
library(stringr)

options(warn = 1)

glue('covidcast MPI runner

Usage:
  {name} [--cpus-per-task=<cores>] [--id-vars=<vars>] --output=<path> <path>
  {name} (-h | --help)
  {name} --version

Options:
  -h --help                 Show this screen.
  --version                 Show version.
  --cpus-per-task=<cores>   How many cores should each task (process) use?
  --output=<path>           Where to save the resulting RDS file
  --id-vars=<vars>          Grouping vars. [default: state]
', name = "MPIrun.R") -> doc

arguments <- docopt(doc, version = 'covidcast MPI runner 0.3')

input_file_path  <- arguments$path
id_vars          <- str_split(arguments[["id-vars"]], ',')[[1]]
output_file_path <- arguments$output
cpus_per_task    <- as.numeric(arguments[['cpus-per-task']])

# Start an MPI cluster, allowing the routine to discover as many tasks
# as are available. This makes it easy to adapt to situations where
# different numbers of tasks are available to use.
startMPIcluster(
  maxcores = cpus_per_task,
  verbose  = TRUE
) -> cl
csize <- clusterSize(cl)
registerDoMPI(cl) # Register the MPI backend with `foreach`'s `%dopar%
cli_alert_success("Initialized MPI cluster")
cli_alert_info("Cluster size is: {.val {csize}}")

# Read input data, make sure you capture the date properly, and be sure that
# there aren't any missing values in the cases/deaths data
read_csv(
  input_file_path,
  col_types = cols(date = col_date(format = "%Y-%m-%d"))
) -> d

assert_that(!any(is.na(d$cases)))
assert_that(!any(is.na(d$deaths)))
assert_that(!any(is.na(d$fracpos)))

cli_alert_success("Read {.file {input_file_path}}, all assertions passed")

d <- local({
  # excluded_states <- c("Puerto Rico")
  included_states <- "Washington"

  # cli_alert_warning("Filtering out the following states:")
  cli_alert_warning("Restricting to the following states:")
  cli_ul(items = included_states)
  # cli_ul(items = excluded_states)

  filter(d, state %in% included_states)
  # filter(d, ! state %in% excluded_states)
})

# This function takes some data for a geographic tract (likely a state), and
# produces a `covidcast` object which, if created without warnings or errors,
# represents a valid model configuration that is ready to be run.
#
# `type_cases`, `type_deaths` should be either "reported" or "occurred"
gen_covidcast_run <- function(d, type_cases = "reported",
                              type_deaths = "reported", N_days_before = 28) {

  d_cases   <- transmute(d, date, observation = cases)
  d_deaths  <- transmute(d, date, observation = deaths)
  d_fracpos <- transmute(d, date, observation = fracpos)
  N_days    <- nrow(d_cases)

  covidcast::covidcast(N_days = N_days,
                       N_days_before = N_days_before) +
    covidcast::input_cases(d_cases,   type = type_cases) +
    covidcast::input_deaths(d_deaths, type = type_deaths) +
    covidcast::input_fracpos(d_fracpos)
}

# Take the data, and split it into smaller tibbles for each state
states_nested <- group_by_at(d, id_vars) %>% nest() %>%
  mutate(config = map(data, gen_covidcast_run))

cli_alert_success("Configuration finished successfully. Summary:")
cli_ul()
pwalk(states_nested, function(data, config, ...) {
  id_vars <- list(...)

  cli_li(paste(id_vars, collapse = '-'))
  ulid <- cli_ul()
  cli_li("Death data: {nrow(data)} obs")
  cli_li("Cases data: {nrow(data)} obs")
  cli_li("Fracpos data: {nrow(data)} obs")
  cli_end(ulid)
})
cli_end()

cli_alert_info("Beginning model executions\n")

foreach(i = seq_along(states_nested$state), # By row
        .inorder = FALSE, # Weird docs suggest this has perf. advantage?
        .combine = 'bind_rows') %dopar% { # Every run returns a 1-row tibble

  cfg                <- states_nested$config[[i]] # config for this run
  # current_state   <- states_nested$state[i] # name of state
  current_group      <- as.list(states_nested[i, id_vars])
  current_group_name <- paste(current_group, collapse = "-")
  safe_run           <- purrr::quietly(covidcast::run) # capture warning msgs
  core               <- Rmpi::mpi.get.processor.name(short = FALSE)
  diagnostic_file    <- glue::glue("logs/stan-{current_group_name}") # Stan logs
  startTime          <- Sys.time()

  cli_alert_info("Started {current_group_name} on {core} at {startTime}")

  safe_run(
    cfg,
    cores = cpus_per_task,
    diagnostic_file = diagnostic_file
  ) -> result
  endTime <- Sys.time()

  cli_alert_info("Finished {current_group_name} on {core} at {endTime}")
  cli_alert_info("Runtime was {prettyunits::pretty_dt(endTime - startTime)}")

  tibble::tibble(!!! current_group, result = list(result))
} -> result

cli_alert_success("Finished model executions\n")

# Serialize to RDS
saveRDS(result, output_file_path)
saveRDS(states_nested, paste0("config_", output_file_path))

cli_alert_success("Saved result to {output_file_path}\n")

# Cleanup the cluster and shut off MPI's communication stuff
closeCluster(cl)
mpi.quit()
