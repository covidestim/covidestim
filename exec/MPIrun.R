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

options(warn = 1)

glue('covidcast MPI runner

Usage:
  {name} [--cpus-per-task=<cores>] --output=<path> <path>
  {name} (-h | --help)
  {name} --version

Options:
  -h --help                 Show this screen.
  --version                 Show version.
  --cpus-per-task=<cores>   How many cores should each task (process) use?
  --output=<path>           Where to save the resulting RDS file
', name = "MPIrun.R") -> doc

arguments <- docopt(doc, version = 'covidcast MPI runner 0.3')

# Start an MPI cluster, allowing the routine to discover as many nodes
# as are available. This makes it easy to adapt to situations where
# different numbers of nodes are available to use.
startMPIcluster(
  maxcores = as.numeric(arguments[['cpus-per-task']]),
  verbose  = TRUE,
  bcast    = FALSE # Should eventually reconsider using this
) -> cl

pname <- mpi.get.processor.name(short = FALSE)
registerDoMPI(cl) # Register the MPI backend with `foreach`'s `%dopar%
sinkWorkerOutput(glue("logs/{pname}.txt")) # Log worker output

cat(glue("Cluster size is: {size}\n", size = clusterSize(cl)))

# Read the states data in, make sure you capture the date properly, and
# be sure that there isn't any missingness in the data
read_csv(
  arguments$path,
  col_types = cols(date = col_date(format = "%Y-%m-%d"))
) -> d

assert_that(!any(is.na(d$incident_cases)))
assert_that(!any(is.na(d$incident_deaths)))

# EXCLUDING ALABAMA, COLORADO, GEORGIA! Each one of these states has
# situations where cumulative deaths or cases go down for at least one day.
d <- filter(d, !state %in% c("Alabama", "Colorado", "Georgia"))

# This function takes some data for one geographic tract, and produces a
# `covidcast` object which, if created without warnings or errors, should
# represent a valid model configuration that is ready to be run.
gen_covidcast_run <- function(d, type = "reported", N_days_before = 10) {

  d_cases  <- transmute(d, date, observation = incident_cases)
  d_deaths <- transmute(d, date, observation = incident_deaths)

  covidcast(N_days = nrow(d_cases),
            N_days_before = N_days_before) +
    input_cases(d_cases,   type = type) +
    input_deaths(d_deaths, type = type)
}

# Take the data, and split it into smaller tibbles for each state
states_nested <- d %>% group_by(state) %>% nest() %>%
  mutate(config = map(data, quietly(gen_covidcast_run)))

cat("Beginning model executions\n")

foreach(i = seq_along(states_nested$state),
        .inorder = FALSE,
        .combine = 'bind_rows') %dopar% {

  cfg             <- states_nested$config[[i]]$result
  current_state   <- states_nested$state[i]
  safe_run        <- purrr::quietly(covidcast::run)
  core            <- Rmpi::mpi.get.processor.name(short = FALSE)
  diagnostic_file <- glue::glue("logs/stan_{current_state}")

  print(glue::glue("[start] {current_state} on {core}"))
  result <- safe_run(cfg, cores = 3, diagnostic_file = diagnostic_file)
  print(glue::glue("[finish] {current_state} on {core}"))

  tibble::tibble(state = current_state, result = list(result))
} -> result

cat("Finished model executions\n")

# Serialize to RDS
saveRDS(result, arguments$output)
saveRDS(states_nested, paste0("config_", arguments$output))

cat(glue("Saved result to {arguments$output}\n"))

# Cleanup the cluster and shut off MPI's communication stuff
closeCluster(cl)
mpi.quit()
