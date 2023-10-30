library(data.table)
library(aod)
library(tidyverse)
library(tictoc)

args = commandArgs(trailingOnly=TRUE)
file_number <- as.integer(args[1])
# file_number = 5 # Remove
 
alt_count_files <- list.files(path="/Genomics/ayroleslab2/scott/git/chromium/data/snakemake_tmp/count_info", pattern = "alternate_counts*", full.names=TRUE)
ref_count_files <- list.files(path="/Genomics/ayroleslab2/scott/git/chromium/data/snakemake_tmp/count_info", pattern = "reference_counts_*", full.names=TRUE)

# fread is faster than read.table
alt_count_data <- read.table(alt_count_files[file_number], header=TRUE, sep="\t", row.names=1)
ref_count_data <- read.table(ref_count_files[file_number], header=TRUE, sep="\t", row.names=1)

both_count_data <- alt_count_data + ref_count_data


info=read_csv('/Genomics/ayroleslab2/scott/git/chromium/metadata/longevity_dna_individual_metadata.csv')
info$name = gsub('2428__','',info$indiv)

info = info[order(info$name),]
both_count_data = both_count_data[,order(colnames(both_count_data))]
alt_count_data = alt_count_data[,order(colnames(alt_count_data))]
colnames(both_count_data) = gsub('X2428__','',colnames(both_count_data))
colnames(alt_count_data) = gsub('X2428__','',colnames(alt_count_data))

# Drop plates
filtered_info <- info %>% filter(!(plate %in% c("Chrom_31_6","Chrom_31_7","Chrom_31_8","Chrom_31_9")))
both_count_data <- both_count_data[,colnames(both_count_data) %in% filtered_info$name]
alt_count_data <- alt_count_data[,colnames(alt_count_data) %in% filtered_info$name]

# If we dropped any individuals, we need to drop them from the info file
info = info[info$name %in% colnames(both_count_data),]

# print shape
print(paste0("Number of individuals: ", dim(both_count_data)[2]))
print(paste0("Number of SNPs: ", dim(both_count_data)[1]))

# Filter snp data to only include SNPs with non-zero alts
filter <- rowSums(alt_count_data,na.rm=TRUE) == 0 | apply(alt_count_data,1,function(x) all(is.na(x)))
alt_count_data <- alt_count_data[!filter,]
both_count_data <- both_count_data[!filter,]

# 
coef_betabin <- matrix(nrow=dim(alt_count_data)[1],ncol=2)

# Number of permutations
B <- 10
plouf <- CJ(it = 1:B) 

number_of_loci = dim(alt_count_data)[1]
# number_of_loci = 100
output_mat = matrix(nrow=number_of_loci, ncol=9)  # Change ncol to 9
colnames(output_mat) <- c("SNP_location", "obs_ctrl","pval_ctrl","obs_crvi","pval_crvi","obs_contrast","pval_contrast","emp_fdr_pvalthresh_05","significant")
tic()
for (i in 1:number_of_loci){
    tryCatch({
        print(paste0("Number " ,i, " of ", number_of_loci, " (", 100*round(i/dim(alt_count_data)[1],5), "%)"))
        # Include SNP location in the output matrix
        snp_location = rownames(both_count_data)[i]

        dt = data.table(depth = unlist(both_count_data[i,],use.names=FALSE), alt_depth = unlist(alt_count_data[i,],use.names=FALSE), sex = info$sex, treatment = info$treatment)
        dt = dt[dt$depth > 0,]

        # Reorder the factor levels to match the order of the contrasts
        dt$treatment <- factor(dt$treatment,levels=c("T0","Control","3mM"))

        # Fit the model -- this is the base model with all three
        mod1 <- tryCatch(betabin(cbind(alt_depth, depth - alt_depth) ~ treatment + sex, ~1,data=dt,method="L-BFGS-B",hessian=TRUE), error=function(x){stop("Model failed")})
        
        # Add try-catch blocks to handle errors in accessing the summary attributes
        obs_ctrl <- tryCatch(attributes(summary(mod1))$Coef[2,1], error=function(x) stop("Model failed"))
        pval_ctrl <- tryCatch(attributes(summary(mod1))$Coef[2,4], error=function(x) stop("Model failed"))

        obs_crvi <- tryCatch(attributes(summary(mod1))$Coef[3,1], error=function(x) stop("Model failed"))
        pval_crvi <- tryCatch(attributes(summary(mod1))$Coef[3,4], error=function(x) stop("Model failed"))

        coef_betabin[i,] <- tryCatch(attributes(summary(mod1))$Coef[2:3,1], error=function(x) stop("Model failed"))

        # Fit the model without the T0 treatment
        dt_no_t0 <- dt[dt$treatment %in% c("Control","3mM")]

        mod2 <- tryCatch(betabin(cbind(alt_depth, depth - alt_depth) ~ treatment + sex, ~1,data=dt_no_t0,method="L-BFGS-B",hessian=TRUE), error=function(x){})
        obs_contrast <- tryCatch(attributes(summary(mod2))$Coef[2,1], error=function(x) stop("Model failed"))
        pval_contrast <- tryCatch(attributes(summary(mod2))$Coef[2,4], error=function(x) stop("Model failed"))
        
        # Setup df for permutations
        f <- function(x) tryCatch(attributes(summary(x))$Coef[2,4], error=function(x) stop("Model failed"))
        dt_no_t0[, .(fitted = f(betabin(cbind(alt_depth, depth - alt_depth) ~ treatment + sex , random = ~1, data = .SD, method="L-BFGS-B",hessian=TRUE)))]

        # Permutations to estimate the FDR
        Ntrial <- plouf[,.(depth = dt_no_t0$depth, alt_depth =dt_no_t0$alt_depth, sex = dt_no_t0$sex, treatment = sample(dt_no_t0$treatment)), by = .(it)]

        res <- Ntrial[, .(tryCatch(f(betabin(cbind(alt_depth, depth - alt_depth) ~ treatment + sex, ~1,data=.SD,method="L-BFGS-B",hessian=TRUE)), error=function(x){ stop("Model failed")})),
            by = .(it)]

        # Calculate the FDR threshold for each locus
        fdr_pval <- tryCatch(quantile(res$V1, 0.1), error=function(x) stop("Model failed"))

        sig <- tryCatch(as.numeric(pval_contrast < fdr_pval), error = function(x) stop("Model failed"))
        
        output_mat[i,] = c(snp_location, obs_ctrl, pval_ctrl, obs_crvi, pval_crvi, obs_contrast, pval_contrast, fdr_pval, sig)
    }, error = function(e) {
        print(paste0("Error at locus ", i))
        print(e)
    })

}
toc()

write.table(output_mat, file=paste0("/Genomics/ayroleslab2/scott/git/chromium/data/20231025_betabin_contrast_w_perm/betabin_perm_output_10perm_fdr010_",file_number,".csv"), sep=",", row.names=FALSE, col.names=TRUE)
