#!/bin/bash
#BATCH --job-name=trim_map_recal_call                   # Name of the job
#SBATCH --output=logs/trim_map_recal_call_%A_%a.out     # STDOUT file
#SBATCH --error=logs/trim_map_recal_call_%A_%a.err      # STDERR file
#SBATCH --nodes=1                                       # Node count
#SBATCH --ntasks=1                                      # Number of tasks across all nodes
#SBATCH --cpus-per-task=4                               # Cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=64G                                       # total memory per node
#SBATCH --array=1-2880                                  # Number of jobs % Max number of jobs to consume
#SBATCH --time=1-00:00:00                               # Run time limit (HH:MM:SS)
#SBATCH --mail-type=all                                 # Email on job start, end, and fault
#SBATCH --mail-user=swwolf@princeton.edu

module load conda
module load samtools

conda activate crvi

export SAMPLEINFO=$(sed "${SLURM_ARRAY_TASK_ID}q;d" list_of_samples.txt)
echo $SAMPLEINFO
sh trim_map_recal_call.sh
