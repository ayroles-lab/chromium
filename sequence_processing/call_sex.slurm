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

for i in {1..2880}
do
    echo $i
    export SLURM_ARRAY_TASK_ID=$i
    export SAMPLEINFO=$(sed "${SLURM_ARRAY_TASK_ID}q;d" list_of_samples.txt)
    export BAM=$(echo ${SAMPLEINFO}_recal.bam)
    export PROCESSORS=$(echo ${SLURM_CPUS_PER_TASK})
    samtools view --threads $PROCESSORS -b /scratch/tmp/swwolf/longevity/$BAM | cut -f 3 | sort | uniq -c | sort -nr | sed -e 's/^ *//;s/ /\t/' | awk 'OFS="\t" {print $2,$1}' | sort -n -k1,1 > /scratch/tmp/swwolf/longevity/${SAMPLEINFO}.chromosome_coverage.tsv &
done

# echo $SAMPLEINFO
# sh trim_map_recal_call.sh
