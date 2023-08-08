snakemake --use-conda --cluster-config cluster_config.json --cluster "sbatch --job-name={cluster.job-name} --output={cluster.output} --error={cluster.error} --nodes={cluster.nodes} --ntasks={cluster.ntasks} --cpus-per-task={cluster.cpus-per-task} --mem={cluster.mem} --time={cluster.time} --mail-type={cluster.mail-type} --mail-user={cluster.mail-user}" --jobs 400 --conda-frontend conda --rerun-incomplete --latency-wait 240


squeue -u $USER | awk '{print $1}' | xargs -n 1 scancel

snakemake --delete-all-output -c24  --rerun-incomplete

rm logs/slurm/*

# Extract the line for each chromosome of interest directly
for chrom in 2L 2R 3L 3R 4 X; dop 
  grep -P "^${chrom}\t" /Genomics/ayroleslab2/scott/git/chromium/data/ref/dmel-all-chromosome-r6.49.fasta.fai
done > chroms.fasta.fai

# Now generate the windows using bedtools
bedtools makewindows -g chroms.fasta.fai -w 350000 > chroms_intervals.bed

for file in *.vcf
do
    tabix -p vcf "${file}.gz" &
done



        r"""
        module load bcftools && \
        mkdir -p {TMP_DIR}/count_info && \
        bcftools query -f '[%CHROM:%POS\\t%SAMPLE\\t%DP\\n]' {input.vcf_file} | sed 's/,/\\t/g' > {TMP_DIR}/count_info/joint_call_{wildcards.chrom}.DP.txt && \
        bcftools query -f '[%CHROM:%POS\\t%SAMPLE\\t%AD\\n]' {input.vcf_file} | sed 's/,/\\t/g' > {TMP_DIR}/count_info/joint_call_{wildcards.chrom}.AD.txt && \
        awk '{{OFS="\\t"; if ($3>0) print $1,$2,$3}}' {TMP_DIR}/count_info/joint_call_{wildcards.chrom}.DP.txt > {output.out_file1} && \
        awk '{{OFS="\\t"; if ($3>0 || $4>0) print $1,$2,$3,$4}}' {TMP_DIR}/count_info/joint_call_{wildcards.chrom}.AD.txt > {output.out_file2} && \
        rm {TMP_DIR}/count_info/joint_call_{wildcards.chrom}.DP.txt && \
        rm {TMP_DIR}/count_info/joint_call_{wildcards.chrom}.AD.txt
        """

bcftools query -f '[%CHROM:%POS\t%SAMPLE\t%DP\n]' /Genomics/ayroleslab2/scott/git/chromium/data/snakemake_results/4.combined.vcf | sed 's/,/\t/g'> /Genomics/ayroleslab2/scott/git/chromium/data/snakemake_tmp/count_info/joint_call_X.DP.txt


plink --vcf /Genomics/ayroleslab2/scott/git/chromium/data/snakemake_results/4.combined.vcf --allow-extra-chr --maf 0.01 --set-missing-var-ids @:#:$1,$2 --indep-pairwise 200 25 0.9 --geno 0.75 --biallelic-only strict --double-id --recode vcf --out /Genomics/ayroleslab2/scott/git/chromium/data/snakemake_results/joint_call_4.SNPs.PASS


        plink --vcf /Genomics/ayroleslab2/scott/git/chromium/data/snakemake_results/4.combined.vcf --allow-extra-chr --maf 0.01 --set-missing-var-ids @:#:$1,$2 --indep-pairwise 200 25 0.9 --geno 0.75 --biallelic-only strict --double-id --recode vcf --out /Genomics/ayroleslab2/scott/git/chromium/data/snakemake_results/joint_call_4.SNPs.PASS                                                                                                                        
        sed -e 's/:/\t/g' /Genomics/ayroleslab2/scott/git/chromium/data/snakemake_results/joint_call_4.SNPs.PASS.prune.in | awk '{OFS="/\t";print $1,$2-1,$2}'> /Genomics/ayroleslab2/scott/git/chromium/data/snakemake_results/joint_call_4.SNPs.PASS.bed                                                                 
                                                                                                                                                               
        bgzip /Genomics/ayroleslab2/scott/git/chromium/data/snakemake_results/joint_call_4.SNPs.PASS.vcf                                                      
        tabix -p vcf /Genomics/ayroleslab2/scott/git/chromium/data/snakemake_results/joint_call_4.SNPs.PASS.vcf.gz                                            
                                                                                                                                                               
        bcftools filter -e "QUAL < 20" --threads 4 -R /Genomics/ayroleslab2/scott/git/chromium/data/snakemake_results/joint_call_4.SNPs.PASS.bed -O z -o /Genomics/ayroleslab2/scott/git/chromium/data/snakemake_results/4.filtered.combined.vcf.gz /Genomics/ayroleslab2/scott/git/chromium/data/snakemake_results/joint_call_4.SNPs.PASS.vcf.gz 


module load R && Rscript --vanilla /Genomics/ayroleslab2/scott/git/chromium/sequence_processing/snakemake/process_counts.R /Genomics/ayroleslab2/scott/git/chromium/data/snakemake_tmp/count_info/joint_call_4.AD.filt.txt /Genomics/ayroleslab2/scott/git/chromium/data/snakemake_tmp/count_info/joint_call_4.DP.filt.txt /Genomics/ayroleslab2/scott/git/chromium/data/snakemake_tmp/count_info/alternate_counts_4.txt /Genomics/ayroleslab2/scott/git/chromium/data/snakemake_tmp/count_info/reference_counts_4.txt