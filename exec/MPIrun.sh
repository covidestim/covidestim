#!/usr/bin/env bash

#SBATCH --cpus-per-task=3
#SBATCH --mail-type=none
#SBATCH --mem-per-cpu=3G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=handsome.dan@yale.edu

module purge
module load miniconda intel/2018b OpenMPI
source activate covidcast

mkdir -p logs

mpirun --bind-to none --display-map Rscript MPIrun.R \
  --id-vars=state,source,smoothing \
  --cpus-per-task="${SLURM_CPUS_PER_TASK-3}" \
  --output "${SLURM_JOB_NAME:-results}.rds" \
  "$1" # The path to the input data

# `--bind-to none` prevents each process from being bound to a core, allowing
#   R's `fork()`s to result in child processes that are bound to different
#   logical cores (since they have been allocated by SLURM for the task)
#
# `--display-map` shows where each worker process is hosted, along with some
#   config information
#
# `--id-vars` specifies which variables in the input data ($1) to use as
#   grouping variables. Each group will be transformed into a model
#   configuration and passed to `covidcast::run()`
