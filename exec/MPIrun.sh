#!/usr/bin/env bash

#SBATCH --cpus-per-task=3
#SBATCH --mail-type=none
#SBATCH --mem-per-cpu=4G

module purge
module load miniconda R OpenMPI
source activate covidcast
mkdir -p logs
                    
echo "Invoking mpirun"

echo mpirun --output-filename logs/mpi \
       --timestamp-output \
  "`dirname $0`/MPIrun.R" \
    --mpi \
    --cpus-per-task=${SLURM_CPUS_PER_TASK:-3} \
    --output "${SLURM_JOB_NAME:-output}.rds" \
    --diag "logs" \
    "$1"

mpirun --output-filename logs/mpi \
       --timestamp-output \
  "`dirname $0`/MPIrun.R" \
    --mpi \
    --cpus-per-task=${SLURM_CPUS_PER_TASK:-3} \
    --output "${SLURM_JOB_NAME:-output}.rds" \
    --diag "logs" \
    "$1"
