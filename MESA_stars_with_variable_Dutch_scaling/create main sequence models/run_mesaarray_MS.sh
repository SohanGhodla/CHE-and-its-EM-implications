#!/bin/bash -e
#SBATCH -J Single_MS
#SBATCH --time=2:00:00
#SBATCH -A uoa03218
#SBATCH --mem-per-cpu=3G
#SBATCH --array=3900-4033
#echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
#SBATCH -o ./Slurm_output/%A_%a.out # STDOUT

srun run_new_Rot_MS.sh "${SLURM_ARRAY_TASK_ID}"

