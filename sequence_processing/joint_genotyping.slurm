#!/bin/bash
#SBATCH --mem=32G
#SBATCH --time=7-00:00:00
#SBATCH --job-name=ImportDB
#SBATCH --cpus-per-task=2
#SBATCH --output="logs/%A_%a.out"
#SBATCH --error="logs/%A_%a.error"
#SBATCH --array=1-128

# run from the inDIR folder where the .list file is located
# requires a mygvcfiles.forImportDB.list with the sample ID of the gvcf files TAB path to the gvcfs

set -e 
date >&2
###########
# joint call by chromosome
###########

cpus=128

ref_DIR=/Genomics/ayroleslab2/scott/git/chromium/data/ref/dmel-all-chromosome-r6.49.fasta
out_DIR=/Genomics/ayroleslab2/scott/git/chromium/data/joint_calls
GATK_DIR=/Genomics/ayroleslab2/scott/git/chromium/tools/gatk
variants_DIR=/Genomics/ayroleslab2/scott/git/chromium/data/ref/dbSNP_Nex_Sep28.19.vcf
# variants_DIR=/Genomics/ayroleslab2/alea/longevity/NovaSeq_Dec18/joint_vcfs/joint_call_n2592.X.PASS.vcf.gz
# $GATK_DIR SplitIntervals -R $ref_DIR --scatter-count 128 -o intervals  

# $GATK_DIR GenomicsDBImport --java-options "-Xmx256G -Xms256G -XX:ParallelGCThreads=32" \
#    --genomicsdb-workspace-path $out_DIR/genomics_dbs/genomicsdb_x \
#    --sample-name-map /Genomics/ayroleslab2/scott/git/chromium/submodules/longevity/vcfs.longevity.map \
#    -L X \
#    --reader-threads 128 \
#    --batch-size 128

# $GATK_DIR CombineGVCFs --java-options "-Xmx256G -XX:ParallelGCThreads=32" \
#    -R $ref_DIR \
#    -V /Genomics/ayroleslab2/scott/git/chromium/submodules/longevity/vcfs.longevity.list \
#    -output $out_DIR/merged_gvcfs_longevity.vcf

# $GATK_DIR GenotypeGVCFs \
#    -R $ref_DIR \
#    -V gendb:///Genomics/ayroleslab2/scott/git/chromium/data/joint_calls/genomic_db_slurm/run1.127 \
#    -O output.vcf.gz --max-alternate-alleles 2

###########
# get filtered SNPs
###########

# https://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set

$GATK_DIR -R $ref_DIR -T SelectVariants --variant $out_DIR/joint_call_X_longevity.vcf -selectType SNP -output $out_DIR/joint_call_X_longevity.SNPs.vcf

rm $out_DIR/joint_call_X_longevity.vcf*

$GATK_DIR -R $ref_DIR -T VariantFiltration -o $out_DIR/joint_call_X_longevity.SNPs.filt.vcf --variant $out_DIR/joint_call_X_longevity.SNPs.vcf --filterExpression "QUAL < 20.0" --filterExpression "QD < 2.0" --filterExpression "FS > 60.0" --filterExpression "MQ < 35.0" --filterExpression "MQRankSum < -12.5" --filterExpression "ReadPosRankSum < -8.0" --filterName "LowQual" --filterName "LowQD" --filterName "HighFS" --filterName "LowMQ" --filterName "LowMQRankSum" --filterName "LowReadPosRankSum"

rm $out_DIR/joint_call_X_longevity.SNPs.vcf*

~/programs/plink_1.90 --vcf $out_DIR/joint_call_X_longevity.SNPs.filt.vcf --allow-extra-chr --maf 0.01 --set-missing-var-ids @:#:\$1,\$2 --indep-pairwise 200 25 0.9 --geno 0.75 --biallelic-only strict --recode vcf --out $out_DIR/joint_call_X_longevity.SNPs.PASS

sed -e s/:/'\t'/g $out_DIR/joint_call_X_longevity.SNPs.PASS.prune.in | awk '{OFS="\t";print $1,$2-1,$2}' > $out_DIR/joint_call_X_longevity.SNPs.PASS.bed

bgzip $out_DIR/joint_call_X_longevity.SNPs.PASS.vcf
tabix -p vcf $out_DIR/joint_call_X_longevity.SNPs.PASS.vcf.gz

bcftools filter -e "QUAL < 20" --threads 4 -R $out_DIR/joint_call_X_longevity.SNPs.PASS.bed -O z -o $out_DIR/joint_call_X_longevity.SNPs.PASS2.vcf.gz $out_DIR/joint_call_X_longevity.SNPs.PASS.vcf.gz

rm $out_DIR/joint_call_X_longevity.SNPs.PASS.*

# keep some unfiltered SNPs and final filtered set
bgzip $out_DIR/joint_call_X_longevity.SNPs.filt.vcf
tabix -p vcf $out_DIR/joint_call_X_longevity.SNPs.filt.vcf.gz

tabix -p vcf $out_DIR/joint_call_X_longevity.SNPs.PASS2.vcf.gz


