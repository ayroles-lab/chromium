library(data.table)
library(aod)
library(tidyverse)
library(rio)
library(ggplot2)

# args = commandArgs(trailingOnly=TRUE)
# file_number =  as.integer(args[1])

list_of_betabin_output_files <- list.files(path="/Genomics/ayroleslab2/scott/git/chromium/data/20231025_betabin_contrast_w_perm/", pattern = "betabin_perm_output_10perm_fdr010*", full.names=TRUE)
betabin_outputs <- lapply(list_of_betabin_output_files, function(x) {import(x, format="csv")})
betabin_outputs <- do.call(rbind, betabin_outputs)
betabin_df <- as.data.frame(betabin_outputs)

p <- ggplot(betabin_df, aes(x=pval_ctrl)) +
        geom_histogram(bins=100) +
        theme_bw() + 
        theme(legend.position="none") + 
        xlab("p-value") + 
        ylab("Number of sites") + 
        ggtitle("Distribution of ctrl p-values from Beta-Binomial model")

ggsave(filename="figures/20231030_6789_removed_betabin_ctrl_pvals_dist.png", plot=p, width=10, height=10, units="in")

betabin_df$chr <- sub(":.*", "", betabin_df$SNP)
betabin_df$BP <- as.numeric(sub(".*:", "", betabin_df$SNP))

drosophila_chromosomes <- c("X", "2R", "2L", "3R", "3L", "4", "Y")

# Select rows with chromosome in Drosophila chromosomes + only keep complete cases
betabin_df <- betabin_df[betabin_df$chr %in% drosophila_chromosomes, ]

# Name columns more descriptively before writing to file
# names(betabin_df) <- c("SNP_ID", "Observed_ctrl", "Control_pvalue", "Observed_crvi", "CRVI_pvalue", "Chromosome", "Basepair")
# write.csv(betabin_df, file="betabin_pval_data.csv", row.names=FALSE)

betabin_df <- betabin_df[complete.cases(betabin_df), ]

# Continue with the script...
betabin_df$pval_ctrl_adj_fdr <- p.adjust(betabin_df$pval_ctrl, method="fdr")
betabin_df$pval_crvi_adj_fdr <- p.adjust(betabin_df$pval_crvi, method="fdr")

gxe_sig <- betabin_df[betabin_df$pval_ctrl_adj_fdr < 0.10 | betabin_df$pval_crvi_adj_fdr < 0.10, ]

# Venn
print("Total number of significant genes across both environments (ctrl/crvi)")
pvals_ctrl_adj_sig <- gxe_sig$SNP[gxe_sig$pval_ctrl_adj_fdr < 0.10]
pvals_crvi_adj_sig <- gxe_sig$SNP[gxe_sig$pval_crvi_adj_fdr < 0.10]

print(length(pvals_ctrl_adj_sig) + length(pvals_crvi_adj_sig))

both_sig <- betabin_df[betabin_df$pval_ctrl_adj_fdr < 0.10 & betabin_df$pval_crvi_adj_fdr < 0.10, ]
ctrl_only_sig <- betabin_df[betabin_df$pval_ctrl_adj_fdr < 0.10 & betabin_df$pval_crvi_adj_fdr >= 0.10, ]
crvi_only_sig <- betabin_df[betabin_df$pval_ctrl_adj_fdr >= 0.10 & betabin_df$pval_crvi_adj_fdr < 0.10, ]

names(both_sig) <- c("SNP", "Observed_ctrl", "Control_pvalue", "Observed_crvi", "CRVI_pvalue", "Chromosome", "Basepair", "Control_FDR_Adj_p-value", "CRVI_FDR_Adj_p-value")
names(ctrl_only_sig) <- c("SNP", "Observed_ctrl", "Control_pvalue", "Observed_crvi", "CRVI_pvalue", "Chromosome", "Basepair", "Control_FDR_Adj_p-value", "CRVI_FDR_Adj_p-value")
names(crvi_only_sig) <- c("SNP", "Observed_ctrl", "Control_pvalue", "Observed_crvi", "CRVI_pvalue", "Chromosome", "Basepair", "Control_FDR_Adj_p-value", "CRVI_FDR_Adj_p-value")

write.csv(both_sig, file="both_environments_sig_pvals.csv", row.names=FALSE)
write.csv(ctrl_only_sig, file="ctrl_only_sig_pvals.csv", row.names=FALSE)
write.csv(crvi_only_sig, file="crvi_only_sig_pvals.csv", row.names=FALSE)


both_sig[1:10, c("SNP", "Control_pvalue", "CRVI_pvalue")]

