#!/usr/bin/env bash

#SBATCH -n 4 -N 4 -t 1:00
#SBATCH --mail-type=none

module purge
module load miniconda

source activate parallel_r

mpirun -np $SLURM_NTASKS \
  ./MPIrun.R \
    --mpi \
    --ncpus-per-task=$SLURM_CPUS_PER_TASK \
    --output "$SLURM_JOB_NAME.rds" \
    "$1"
