#!/bin/bash
#SBATCH --mem=900G
#SBATCH --time=7-00:00:00
#SBATCH --job-name=ImportDB
#SBATCH --cpus-per-task=4
#SBATCH --output="logs/%A_%a.out"
#SBATCH --error="logs/%A_%a.error"
#SBATCH --array=1-128
# run from the inDIR folder where the .list file is located
# requires a mygvcfiles.forImportDB.list with the sample ID of the gvcf files TAB path to the gvcfs

set -e 
date >&2

# refgenome=/Genomics/grid/users/lamaya/genomes/dmel_genome/dmel-all-chromosome-r6.14.fa
# vcfout=bigproject_batch1.vcf
# find ~+ -type f -name . > mapped_vcfs.list 

logfile=log.importDB.longevity
inDIR=/scratch/tmp/lamaya/novaseq_gDNA
gvcfs=/Genomics/ayroleslab2/scott/git/chromium/submodules/longevity/vcfs.longevity.map
dbpath=/Genomics/ayroleslab2/scott/git/chromium/data/joint_calls/genomic_db_slurm/filtered/
chrlist=/Genomics/ayroleslab2/scott/git/chromium/submodules/longevity/interval_files.list
gatk=/Genomics/ayroleslab2/scott/git/chromium/tools/gatk


# ref_DIR=/Genomics/ayroleslab2/scott/git/chromium/data/ref/dmel-all-chromosome-r6.49.fasta
# out_DIR=/Genomics/ayroleslab2/scott/git/chromium/data/joint_calls
# gatk=/Genomics/ayroleslab2/scott/git/chromium/tools/gatk
# variants_DIR=/Genomics/ayroleslab2/scott/git/chromium/data/ref/dbSNP_Nex_Sep28.19.vcf

source /Genomics/ayroleslab2/scott/git/chromium/submodules/longevity/dbimport.sh ${gvcfs} ${logfile} ${inDIR} ${dbpath} ${chrlist} ${gatk}


date >&2