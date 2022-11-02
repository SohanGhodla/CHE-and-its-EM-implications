#!/bin/bash -e
#SBATCH -J MESA_models
#SBATCH --time=10:00:00
#SBATCH -A uoa03218
#SBATCH --mem-per-cpu=3G
#SBATCH --array=0-406
#echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
#SBATCH -o ./Slurm_output/%A_%a.out # STDOUT

srun run_Rot.sh "${SLURM_ARRAY_TASK_ID}"

