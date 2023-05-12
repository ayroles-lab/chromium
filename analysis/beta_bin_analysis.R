library(data.table)
library(aod)
library(tidyverse)
library(rio)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
file_number =  <- as.integer(args[1])

list_of_betabin_output_files <- list.files(path="/Genomics/ayroleslab2/scott/git/chromium/data", pattern = ".*betabin_gc_.*", full.names=TRUE)
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

ggsave(filename="figures/betabin_ctrl_pvals_dist.png", plot=p, width=10, height=10, units="in")

betabin_df$chr <- sub(":.*", "", betabin_df$SNP)
betabin_df$BP <- as.numeric(sub(".*:", "", betabin_df$SNP))

drosophila_chromosomes <- c("X", "2R", "2L", "3R", "3L", "4", "Y")

# Select rows with chromosome in Drosophila chromosomes + only keep complete cases
betabin_df <- betabin_df[betabin_df$chr %in% drosophila_chromosomes, ]
betabin_df <- betabin_df[complete.cases(betabin_df), ]


betabin_df$pval_ctrl_adj_fdr <- p.adjust(betabin_df$pval_ctrl, method="fdr")
betabin_df$pval_crvi_adj_fdr <- p.adjust(betabin_df$pval_crvi, method="fdr")

gxe_sig <- betabin_df[betabin_df$pval_ctrl_adj_fdr < 0.05 | betabin_df$pval_crvi_adj_fdr < 0.05, ]
dim(gxe_sig)

print("Total number of significant genes across both environments (ctrl/crvi)"
print(length(pvals_ctrl_adj_sig) + length(pvals_crvi_adj_sig))

write.csv(betabin_df, file="betabin_df.csv")
betabin_df$logp <- -log10(betabin_df$pval_ctrl)


alpha <- 0.05 / length(betabin_df$pval_ctrl)
threshold <- -log10(0.05)

library(qqman)

png("figures/ctrl_pvals_qqplot.jpg")
qq(betabin_df$pval_ctrl)
dev.off()

png("figures/ctrl_pvals_qqplot.jpg")
qq(betabin_df$pval_crvi)
dev.off()
