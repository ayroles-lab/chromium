#!/bin/sh

# NOTE: "${SAMPLEINFO}" was replaced with the sample name

#######
# set parameters and directories
#######




cpus=4
ref_DIR=/Genomics/ayroleslab2/scott/git/chromium/data/ref/dmel-all-chromosome-r6.49.fasta
fastq_DIR=/Genomics/ayroleslab2/scott/git/chromium/data/raw
picard_DIR=/Genomics/ayroleslab2/scott/git/chromium/tools/picard.jar
sambamba_DIR=sambamba
GATK_DIR=/Genomics/ayroleslab2/scott/git/chromium/tools/gatk
variants_DIR=/Genomics/ayroleslab2/scott/git/chromium/data/ref/dbSNP_Nex_Sep28.19.vcf

tmp_DIR=/scratch/tmp/swwolf/longevity

# ${SAMPLEINFO}-read-1.fastq.gz
# ${SAMPLEINFO}-read-4.fastq.gz

R1=$fastq_DIR/${SAMPLEINFO}-read-1.fastq.gz
R2=$fastq_DIR/${SAMPLEINFO}-read-4.fastq.gz

R1_trim=$tmp_DIR/${SAMPLEINFO}.trim.R1.fastq.gz
R2_trim=$tmp_DIR/${SAMPLEINFO}.trim.R2.fastq.gz

#######
# trim
#######
# TODO: Double check this!
cutadapt -e 0.1 --overlap 2 -a AGATCGGAAGAG -A AGATCGGAAGAG --minimum-length=20 --trim-n -j 0 -o $R1_trim -p $R2_trim $R1 $R2

#######
# map
#######
bwa mem -M -t $cpus $ref_DIR $R1_trim $R2_trim | samtools view -Sbq 1 | samtools sort -@ $cpus -m 3G > $tmp_DIR/${SAMPLEINFO}_sort_uniq.bam
samtools index $tmp_DIR/${SAMPLEINFO}_sort_uniq.bam
samtools view $tmp_DIR/${SAMPLEINFO}_sort_uniq.bam -b 2L 2R 3L 3R X 4 > $tmp_DIR/${SAMPLEINFO}_sort_uniq_filt_to_chr.bam
samtools sort -@ $cpus -m 3G $tmp_DIR/${SAMPLEINFO}_sort_uniq_filt_to_chr.bam > $tmp_DIR/${SAMPLEINFO}_sort_uniq_filt_to_chr_sort.bam
samtools index $tmp_DIR/${SAMPLEINFO}_sort_uniq_filt_to_chr_sort.bam
# samtools idxstats $tmp_DIR/${SAMPLEINFO}_sort_uniq_filt_to_chr_sort.bam | cut -f 1,3
# samtools view -b $tmp_DIR/${SAMPLEINFO}_sort_uniq_filt_to_chr_sort.bam

rm $R1_trim
rm $R2_trim

#######
# dedup
#######
java -Xmx40g -jar $picard_DIR MarkDuplicates I=$tmp_DIR/${SAMPLEINFO}_sort_uniq_filt_to_chr_sort.bam O=$tmp_DIR/${SAMPLEINFO}_sort_uniq_dedup.bam M=$tmp_DIR/${SAMPLEINFO}_sort_uniq_dedup.txt

$sambamba_DIR index -t $cpus $tmp_DIR/${SAMPLEINFO}_sort_uniq_dedup.bam

#######
# add read groups
#######
java -jar $picard_DIR AddOrReplaceReadGroups I=$tmp_DIR/${SAMPLEINFO}_sort_uniq_dedup.bam O=$tmp_DIR/${SAMPLEINFO}_sort_uniq_dedup_rg.bam SO=coordinate RGLB=${SAMPLEINFO} RGPL=illumina RGPU=longevity RGSM=${SAMPLEINFO} 

$sambamba_DIR index -t $cpus $tmp_DIR/${SAMPLEINFO}_sort_uniq_dedup_rg.bam

#######
# base recalibration and application
#######
# java -Xmx40g -jar $GATK_DIR BaseRecalibrator -R $ref_DIR -I $tmp_DIR/${SAMPLEINFO}_sort_uniq_dedup_rg.bam -output $tmp_DIR/${SAMPLEINFO}.recal.table -known-sites $variants_DIR
$GATK_DIR BaseRecalibrator -R $ref_DIR -I $tmp_DIR/${SAMPLEINFO}_sort_uniq_dedup_rg.bam -output $tmp_DIR/${SAMPLEINFO}.recal.table -known-sites $variants_DIR
# $GATK_DIR PrintReads -R $ref_DIR -I $tmp_DIR/${SAMPLEINFO}_sort_uniq_dedup_rg.bam -BQSR $tmp_DIR/${SAMPLEINFO}.recal.table -o $tmp_DIR/${SAMPLEINFO}_recal.bam
$GATK_DIR ApplyBQSR -R $ref_DIR -I $tmp_DIR/${SAMPLEINFO}_sort_uniq_dedup_rg.bam --bqsr-recal-file $tmp_DIR/${SAMPLEINFO}.recal.table -O $tmp_DIR/${SAMPLEINFO}_recal.bam
$sambamba_DIR index -t $cpus $tmp_DIR/${SAMPLEINFO}_recal.bam

#######
# individual variant calling
#######
$GATK_DIR HaplotypeCaller -R $ref_DIR -I $tmp_DIR/${SAMPLEINFO}_recal.bam -ERC GVCF -output $tmp_DIR/${SAMPLEINFO}.recal.vcf  -max-alternate-alleles 2 --native-pair-hmm-threads $cpus

# $GATK_DIR HaplotypeCaller -R $ref_DIR -I $tmp_DIR/${SAMPLEINFO}_sort_uniq_dedup_rg.bam -ERC GVCF -output $tmp_DIR/${SAMPLEINFO}_sort_uniq_dedup_rg.vcf -max-alternate-alleles 2 --native-pair-hmm-threads 128

bgzip $tmp_DIR/${SAMPLEINFO}.recal.vcf
tabix -p vcf $tmp_DIR/${SAMPLEINFO}.recal.vcf.gz

#######
# clean up
#######
rm $tmp_DIR/${SAMPLEINFO}_sort_uniq*


