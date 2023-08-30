#!/bin/bash
#SBATCH --mem=256G
#SBATCH --time=7-00:00:00
#SBATCH --job-name=betabin_perm_20_fdr05
#SBATCH --cpus-per-task=1
#SBATCH --output="logs/beta_binomial_perm_%A_%a.out"
#SBATCH --error="logs/beta_binomial_perm_%A_%a.error"
#SBATCH --array=1-6

module load R

Rscript beta_binomial_perm.R $SLURM_ARRAY_TASK_ID