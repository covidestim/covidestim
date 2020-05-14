#!/usr/bin/env Rscript
library(Rmpi)
suppressPackageStartupMessages(library(doMPI))
library(docopt)
library(readr)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(covidcast)
library(assertthat)
library(purrr, warn.conflicts = FALSE)
library(glue,  warn.conflicts = FALSE)

glue('covidcast MPI runner

Usage:
  {name} [--mpi] [--cpus-per-task=<cores>] --output=<path> <path>
  {name} (-h | --help)
  {name} --version

Options:
  -h --help                 Show this screen.
  --version                 Show version.
  --mpi                     Attempt to run in parallel using MPI
  --cpus-per-task=<cores>   How many cores should each task (process) use?
  --output=<path>           Where to save the resulting RDS file
', name = "MPIrun.R") -> doc

arguments <- docopt(doc, version = 'covidcast MPI runner 0.2')
print(arguments)

if (arguments$mpi) {
  # Start an MPI cluster, allowing the routine to discover as many nodes
  # as are available. This makes it easy to adapt to situations where
  # different numbers of nodes are available to use.
  cl <- startMPIcluster()
  registerDoMPI(cl) # Register the MPI backend with `foreach`'s `%dopar%
}

# Read the states data in, make sure you capture the date properly, and
# be sure that there isn't any missingness in the data
read_csv(
  arguments$path,
  col_types = cols(
    date = col_date(format = "%Y-%m-%d")
  )
) -> d

assert_that(!any(is.na(d$incident_cases)))
assert_that(!any(is.na(d$incident_deaths)))

######## EXCLUDING ALABAMA, COLORADO, GEORGIA! Each one of these states
######## has situations where cumulative deaths or cases go down for at
######## least one day.
d <- filter(d, !state %in% c("Alabama", "Colorado", "Georgia"))
########
########
########

# This function takes some data for one geographic tract, and produces
# a `covidcast` object which, if created without warnings or errors,
# should represent a valid model configuration that is ready to be run.
gen_covidcast_run <- function(d, type = "reported", N_days_before = 10) {

  d_cases  <- transmute(d, date, observation = incident_cases)
  d_deaths <- transmute(d, date, observation = incident_deaths)

  assert_that(identical(nrow(d_cases), nrow(d_deaths)))

  covidcast(N_days = nrow(d_cases),
            N_days_before = N_days_before) +
    input_cases(d_cases,   type = type) +
    input_deaths(d_deaths, type = type)
}

# Take the data, and split it into smaller tibbles for each state
d %>% group_by(state) %>% nest() -> states_nested

mutate(
  # Right now, running on just 3 states becasue we're not on the cluster
  states_nested[1:3,],

  # Create the `covidconfig` object for each state
  config = map(data, quietly(gen_covidcast_run)),

  # Execute each object by calling `run` on it, through MPI
  result = foreach(
    cfg_item = config, current_state = state
  ) %dopar% {
    # core <- mpi.get.processor.name(short = FALSE)
    run(cfg_item$result, cores = arguments[['cpus-per-tas']])
  }
) -> result

# Pretty print the `result` tibble
print(result)

# Serialize to RDS
saveRDS(result, arguments$output)

if (arguments$mpi) {
  # Cleanup the cluster and shut off MPI's communication stuff
  closeCluster(cl)
  mpi.quit()
}
