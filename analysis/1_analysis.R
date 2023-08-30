library(data.table)
library(tidyverse)
library(rio)
library(ggplot2)


list_of_betabin_output_files <- list.files(path="/Genomics/ayroleslab2/scott/git/chromium/data/20230823_betabin_contrast_w_perm/", pattern = "betabin_perm_output_10perm_fdr010*", full.names=TRUE)
betabin_outputs <- lapply(list_of_betabin_output_files, function(x) {import(x, format="csv")})
betabin_outputs <- do.call(rbind, betabin_outputs)
betabin_df <- as.data.frame(betabin_outputs, stringsAsFactors = FALSE)

# 1. FDR correction on the p-values
betabin_df$pval_ctrl_adj <- p.adjust(betabin_df$pval_ctrl, method = "fdr")
betabin_df$pval_contrast_adj <- p.adjust(betabin_df$pval_contrast, method = "fdr")

# 2. Identify values that are not significant in pval_ctrl but are significant in pval_contrast.
# Let's use a significance threshold of 0.05 (adjust as needed).
threshold <- 0.05

significant_contrast_df <- betabin_df[
  (betabin_df$pval_ctrl_adj > threshold) & (betabin_df$pval_contrast_adj <= threshold),
]

# 70720 SNPs...
dim(significant_contrast_df)
