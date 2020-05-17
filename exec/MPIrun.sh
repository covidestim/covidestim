#!/usr/bin/env bash

#SBATCH --cpus-per-task=3
#SBATCH --mail-type=none
#SBATCH --mem-per-cpu=4G
#SBATCH --mail-type=ALL

module purge
module load miniconda
source activate covidcast

mkdir -p logs
                    
echo "Invoking mpirun"

mpirun -v --output-filename logs/mpi --timestamp-output \
  Rscript MPIrun.R \
    --cpus-per-task=${SLURM_CPUS_PER_TASK:-3} \
    --output "${SLURM_JOB_NAME:-output}.rds" \
    --diag "logs/" \
    "$1"
