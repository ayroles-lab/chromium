#!/bin/bash
#SBATCH --mem=64G
#SBATCH --time=3-00:00:00
#SBATCH --job-name=betabin_perm_10k
#SBATCH --cpus-per-task=2
#SBATCH --output="logs/beta_binomial_perm_%A_%a.out"
#SBATCH --error="logs/beta_binomial_perm_%A_%a.error"
#SBATCH --array=1-128

module load R

Rscript beta_binomial_perm.R $SLURM_ARRAY_TASK_ID