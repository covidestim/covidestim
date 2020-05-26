#!/usr/bin/env bash

#SBATCH -p covid
#SBATCH -A covid
#SBATCH -J run4-allstates-smoothed
#SBATCH --cpus-per-task=3
#SBATCH --mail-type=none
#SBATCH --mem-per-cpu=1G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marcus.russi@yale.edu

module purge
module load miniconda
source activate covidcast

mkdir -p logs
                    
echo "Invoking mpirun"

# mpirun -v --output-filename logs/mpi --timestamp-output \
mpirun -v --output-filename logs/mpi \
  Rscript MPIrun.R \
    --cpus-per-task=${SLURM_CPUS_PER_TASK:-3} \
    --output "${SLURM_JOB_NAME:-output}.rds" \
    "$1" # The path to the input data
