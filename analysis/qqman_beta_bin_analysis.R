library(data.table)
library(aod)
library(tidyverse)
library(rio)
library(ggplot2)
library(qqman)
library(ggman)

args = commandArgs(trailingOnly=TRUE)
file_number <- as.integer(args[1])
file_number = 5 # Remove

list_of_betabin_output_files <- list.files(path="/Genomics/ayroleslab2/scott/git/chromium/data/betabin_contrast_w_perm/", pattern = "betabin_perm_output_10perm_fdr010*", full.names=TRUE)
betabin_outputs <- lapply(list_of_betabin_output_files, function(x) {import(x, format="csv")})
betabin_outputs <- do.call(rbind, betabin_outputs)
betabin_df <- as.data.frame(betabin_outputs, stringsAsFactors = FALSE)

betabin_df$chr <- sub(":.*", "", betabin_df$SNP)
betabin_df$BP <- as.numeric(sub(".*:", "", betabin_df$SNP))

drosophila_chromosomes <- c("X", "2R", "2L", "3R", "3L", "4", "Y")

# Select rows with chromosome in Drosophila chromosomes + only keep complete cases
betabin_df <- betabin_df[betabin_df$chr %in% drosophila_chromosomes, ]

pval_cols <- c("pval_ctrl", "pval_crvi")
window_sizes <- c(0, 150, 500)

for (window_size in window_sizes) {
  if (window_size > 0) {
    betabin_df_windowed <- betabin_df %>%
      mutate(window = cut(BP, breaks = seq(min(BP), max(BP), by = window_size))) %>%
      group_by(chr, window) %>%
      filter(pval_ctrl == min(pval_ctrl), pval_crvi == min(pval_crvi)) %>%
      ungroup()
  } else {
    betabin_df_windowed <- betabin_df
  }

  betabin_df_windowed$chr <- as.numeric(factor(betabin_df_windowed$chr, levels=drosophila_chromosomes))
  for (pval_col in pval_cols) {
    betabin_df_windowed[[pval_col]] <- sapply(betabin_df_windowed[[pval_col]], function(x) ifelse(is.finite(x), x, NA))

    # Manhattan plot
    ggman(betabin_df_windowed, chrom="chr", bp="BP", pvalue="obs_crvi", snp="SNP_location", col=c("blue4", "orange3"))
    ggsave(filename=paste0("figures/manhattan_plot_", pval_col, "_window_", window_size, ".jpg"), width=10, height=10, units="in")

    # QQ plot
    png(paste0("figures/", pval_col, "_qqplot_window_", window_size, ".jpg"))
    qq(betabin_df_windowed[[pval_col]])
    dev.off()
  }
}


