#########################
####MiSeq QC pipeline####
#########################

#General notes
#Note that there are a few steps here (cd'ing around, moving outfiles, etc.) that are omitted, so this should not run as one code chunk
#However, it should serve as a guide to following my QC pipeline

#The basic structure is that in my chromium miseq folder (/scratch/tmp/ed7982/chromium/miseq), there are a series of ordered folders with numeric prefixes (one per step).
#These are 00_resources  01_raw_fastqs  02_demultiplexed_fastqs  03_trimmed_galore  04_raw_bam  05_cleaned_sorted_bam  06_deduped_bam  07_samtools_stats
#00_resources should contain the genome reference and a one-column file of sample names (sample_ids.txt); everything else is part of the pipeline

####Demultiplexing####

#Get files from HTSEQ
wget -r -nH --cut-dirs=2 https://htseq.princeton.edu/tmp/TnBqsapOSuw9gLjCH/

#Make tsv with sample names in first column and two barcodes joined by a hyphen (i7-i5) in second, for demultiplexing
awk '{print $1 "\t" $2 "-" $3}' < longevity_barcodes_all.tsv > barcodes_fmtd.fil

#Get fastq-multx from github
git clone https://github.com/brwnj/fastq-multx
cd fastq-multx
make

#Run the demultiplexing (on files in 01_raw_fastqs)
/Genomics/grid/users/ed7982/fastq-multx/fastq-multx -B barcodes_fmtd.fil \
-m 1 \
4588__VA_chrom_longev-for-318-cycles-000000000-DHCH6_1_Read_2_Index_Read_passed_filter.fastq.gz \
4588__VA_chrom_longev-for-318-cycles-000000000-DHCH6_1_Read_3_Index_Read_passed_filter.fastq.gz \
4588__VA_chrom_longev-for-318-cycles-000000000-DHCH6_1_Read_1_passed_filter.fastq.gz \
4588__VA_chrom_longev-for-318-cycles-000000000-DHCH6_1_Read_4_passed_filter.fastq.gz \
-o n/a \
-o n/a \
-o /scratch/tmp/ed7982/chromium/miseq/02_demultiplexed_fastqs/%_R1.fastq \
-o /scratch/tmp/ed7982/chromium/miseq/02_demultiplexed_fastqs/%_R2.fastq


####Trim Galore (https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)####

#Trim Galore cleans up reads; it's a handy wrapper for fastqc and cutadapt

#Start by getting a file of sample IDs from the raw fastqs folder
awk '{print $1}' < barcodes_fmtd.fil > ../00_resources/sample_ids.txt

#Run Trim Galore (this code chunk should run in 03_trimmed_galore)
while IFS='$\t' read -r id
do
sbatch --wrap="\
/Genomics/grid/users/ed7982/TrimGalore-0.6.6/trim_galore \
--path_to_cutadapt /Genomics/grid/users/ed7982/.local/bin/cutadapt \
--paired \
--basename ${id} \
/scratch/tmp/ed7982/chromium/miseq/02_demultiplexed_fastqs/${id}_R1.fastq \
/scratch/tmp/ed7982/chromium/miseq/02_demultiplexed_fastqs/${id}_R2.fastq"
done < ../00_resources/sample_ids.txt

####BWA mem for alignment####

#Make sure you've downloaded the latest genome reference from Flybase to 00_resources (not shown) and indexed it as below
bwa index dmel-all-chromosome-r6.47.fasta

#Run BWA to align
while IFS='$\t' read -r sample
do
sbatch --wrap="\
module load samtools
fq_fileR1=\"/scratch/tmp/ed7982/chromium/miseq/03_trimmed_galore/${sample}_val_1.fq\"
fq_fileR2=\"/scratch/tmp/ed7982/chromium/miseq/03_trimmed_galore/${sample}_val_2.fq\"
bam_file=\"/scratch/tmp/ed7982/chromium/miseq/04_raw_bam/${sample}.bam\"
bwa_db=\"/scratch/tmp/ed7982/chromium/miseq/00_resources/dmel-all-chromosome-r6.47.fasta\"
bwa mem -M -t 4 -R \"@RG\tID:${sample}\tLB:${sample}\tPL:Illumina\tPU:${sample}\tSM:${sample}\" \$bwa_db \$fq_fileR1 \$fq_fileR2 | samtools view -b > \$bam_file"
done < /scratch/tmp/ed7982/chromium/miseq/00_resources/sample_ids.txt

####Picard for cleaning, sorting, deduping####

#Make sure you have a working path to Picard

while IFS='$\t' read -r sample
do
sbatch --wrap="\
java -jar /Genomics/grid/users/ed7982/picard.jar CleanSam \
    I=/scratch/tmp/ed7982/chromium/miseq/04_raw_bam/${sample}.bam \
    O=/dev/stdout \
    QUIET=true \
    COMPRESSION_LEVEL=0 |
java -jar /Genomics/grid/users/ed7982/picard.jar SortSam \
    I=/dev/stdin \
    O=/scratch/tmp/ed7982/chromium/miseq/05_cleaned_sorted_bam/${sample}.clean.sort.bam \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true 
java -jar /Genomics/grid/users/ed7982/picard.jar MarkDuplicates \
    I=/scratch/tmp/ed7982/chromium/miseq/05_cleaned_sorted_bam/${sample}.clean.sort.bam \
    O=/scratch/tmp/ed7982/chromium/miseq/06_deduped_bam/${sample}.clean.sort.dedup.bam \
    ASSUME_SORTED=true \
    CREATE_INDEX=true \
    M=/scratch/tmp/ed7982/chromium/miseq/06_deduped_bam/${sample}_dedup_metrics.txt"
done < /scratch/tmp/ed7982/chromium/miseq/00_resources/sample_ids.txt


####samtools stats####

#We run samtools stats on our final bams, and then extract the most important metrics from the file (the SN section from http://www.htslib.org/doc/samtools-stats.html)

module load samtools
while IFS='$\t' read -r sample
do
echo "${sample}" $(samtools stats /scratch/tmp/ed7982/chromium/miseq/06_deduped_bam/${sample}.clean.sort.dedup.bam | grep ^SN | cut -f 2- | sed 's/^.*://' | sed 's/#.*$//') > /scratch/tmp/ed7982/chromium/miseq/07_samtools_stats/temp_string.txt
cat /scratch/tmp/ed7982/chromium/miseq/07_samtools_stats/temp_string.txt >> /scratch/tmp/ed7982/chromium/miseq/07_samtools_stats/samtools_stats.txt
done < /scratch/tmp/ed7982/chromium/miseq/00_resources/sample_ids.txt
rm /scratch/tmp/ed7982/chromium/miseq/07_samtools_stats/temp_string.txt

####Adding headers####
#We need to manually add in the headers to the samtools stats text document. I do this in R but you can do it in anything
#After downloading samtools_stats.txt, I run the following in R on my local machine
library(tidyverse)
samtools_stats <- read.table("~/Downloads/samtools_stats.txt", quote="\"", comment.char="")
colnames(samtools_stats) <- c("sample", "raw_total_sequences", "filtered_sequences", "sequences", "is_sorted", "1st_fragments","last_fragments","reads_mapped","reads_mapped_and_paired","reads_unmapped","reads_properly_paired","reads_paired","reads_duplicated","reads_MQ0","reads_QC_failed","non_primary_alignments","total_length","total_first_fragment_length","total_last_fragment_length","bases_mapped","bases_mapped_cigar","bases_trimmed","bases_duplicated","mismatches","error_rate","average_length","average_first_fragment_length","average_last_fragment_length","maximum_length","maximum_first_fragment_length","maximum_last_fragment_length","average_quality","insert_size_average","insert_size_standard_deviation","inward_oriented_pairs","outward_oriented_pairs","pairs_with_other_orientation","pairs_on_different_chromosomes","percentage_of_properly_paired_reads")
write_tsv(samtools_stats, "~/Downloads/chromium_stats_with_colnames.txt")

#From here, plot and QC away!