#!/bin/bash Rscript

# the full list of sites were broken up into parts and run in parallel
# for f in `seq 1 5000 280000`; do cat bb.sh | sed -e s/STARTNUMBER/$f/g > bb.$f.sh; done
# for f in `seq 1 5000 280000`; do cat bb2.sh | sed -e s/STARTNUMBER/$f/g > bb2.$f.sh; done
# rm commands1.sh; touch commands1.sh; for f in `seq 1 5000 280000`; do echo "sh bb2.$f.sh" >> commands1.sh; done


args = commandArgs(trailingOnly=TRUE)
file_number =  <- as.integer(args[1])

library(data.table)
library(aod)
library(tidyverse)
library(rio)


# list_of_glm_output_files <- list.files(path="/Genomics/ayroleslab2/scott/git/chromium/submodules/longevity", pattern = ".*_glm_permutation_pvals.*", full.names=TRUE)
# glm_outputs <- import(list_of_glm_output_files, format="csv")
list_of_glm_output_files <- list.files(path="/Genomics/ayroleslab2/scott/git/chromium/submodules/longevity", pattern = ".*_glm_permutation_pvals.*", full.names=TRUE)
glm_outputs <- lapply(list_of_glm_output_files, function(x) {import(x, format="csv")})
glm_outputs <- do.call(rbind, glm_outputs)
glm_df <- as.data.frame(glm_outputs)

library(ggplot2)

p <- ggplot(glm_df, aes(x=pvalues_ctrl)) +
		geom_histogram(bins=100) +
		theme_bw() + 
		theme(legend.position="none") + 
		xlab("p-value") + 
		ylab("Number of sites") + 
		ggtitle("Distribution of p-values from GLM permutation test")

ggsave(filename="tmp.png", plot=p, width=10, height=10, units="in")

glm_df$chr <- sub(":.*", "", glm_df$snp)
glm_df$BP <- as.numeric(sub(".*:", "", glm_df$snp))

drosophila_chromosomes <- c("X", "2R", "2L", "3R", "3L", "4", "Y")

# Select rows with chromosome information in Drosophila chromosomes
glm_df <- glm_df[glm_df$chr %in% drosophila_chromosomes, ]
glm_df$logp <- -log10(glm_df$pvalues_ctrl)

library(ggman)
install.packages("ggman")

alpha <- 0.05 / length(glm_df$pvalues_ctrl)
threshold <- -log10(0.05)
glm_df$pvalues_ctrl[glm_df$pvalues_ctrl == 0 ] <- 1e-10

library(qqman)
glm_df$chr <- as.numeric(as.factor(glm_df$chr))
png("tmp.png")
# manhattan(glm_df, chr="chr", snp="pos", p="pvalues_ctrl", main="QQ plot of p-values from GLM permutation test", xlab="Expected -log10(p-value)", ylab="Observed -log10(p-value)")
qq(glm_df$pvalues_ctrl)
dev.off()

# glm_df |> group_by(pvalues_ctrl) |> summarise(n=n())
# ggsave(filename="tmp.png", plot=p, width=10, height=10, units="in")


# plot a qq plot with ggplot