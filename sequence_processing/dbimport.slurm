#!/bin/bash
#SBATCH --mem=900G
#SBATCH --time=7-00:00:00
#SBATCH --job-name=ImportDB
#SBATCH --cpus-per-task=24
#SBATCH --output="logs/dbimport_%A_%a.out"
#SBATCH --error="logs/dbimport_%A_%a.err"
#SBATCH --array=1-6
# run from the inDIR folder where the .list file is located
# requires a mygvcfiles.forImportDB.list with the sample ID of the gvcf files TAB path to the gvcfs

set -e 
date >&2

# refgenome=/Genomics/grid/users/lamaya/genomes/dmel_genome/dmel-all-chromosome-r6.14.fa
# find ~+ -type f -name . > mapped_vcfs.list 

logfile=log.importDB.crvi
inDIR=none
gvcfs=/Genomics/ayroleslab2/scott/git/chromium/sequence_processing/vcfs.longevity.map
dbpath=/Genomics/ayroleslab2/scott/git/chromium/data/joint_calls/genomic_db_slurm/20230519_run/
chrlist=/Genomics/ayroleslab2/scott/git/chromium/sequence_processing/chromosomes.list
gatk=/Genomics/ayroleslab2/scott/git/chromium/tools/gatk

source dbimport.sh ${gvcfs} ${logfile} ${inDIR} ${dbpath} ${chrlist} ${gatk}

date >&2