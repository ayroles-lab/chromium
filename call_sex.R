library(tidyverse)
library(ggplot2)


indiv_chr_counts = grep("*chromosome_coverage", list.files("/scratch/tmp/swwolf/longevity", full.names = TRUE), value = TRUE)

df = data.frame()

for (i in 1:length(indiv_chr_counts)) {
  chr_counts = read_tsv(indiv_chr_counts[i], col_names=FALSE) %>% 
    rename(chr = X1, count = X2)
    freq_x = chr_counts$count[chr_counts$chr == "X"]/sum(chr_counts$count[chr_counts$chr != "X"])
    name = gsub(indiv_chr_counts[i], pattern = "*.chromosome_coverage.tsv", replacement = "")
    name = gsub(name, pattern = "/scratch/tmp/swwolf/longevity/", replacement = "")
    df = rbind(df, data.frame(indiv = name, freq_x = freq_x))
}


threshold = 0.175

x_read_freq_dist = ggplot(df, aes(x = freq_x)) +
    geom_vline(xintercept = threshold, linetype = "dashed", color = "red") +
    geom_histogram(bins = 100) +
    theme_minimal() + theme(legend.position = "none") +
    xlab("Frequency of X chromosome reads") + ylab("Number of individuals") +
    ggtitle("Distribution of X chromosome read frequencies")
ggsave("figures/x_read_freq_dist.png", x_read_freq_dist, width = 16, height = 9, units = "in", dpi = 300)


df[df$freq_x < threshold,"sex"]= "male"
df[df$freq_x >= threshold,"sex"] = "female"

library(stringr)

pattern <- "(?<=__)Chrom_\\d+_\\d+"

df$plate = str_extract(df$indiv, pattern)
plate_map = read_csv("data/plate_treatment_mapping.csv", col_names = TRUE)

library(plyr)

df$treatment <- mapvalues(df$plate, plate_map$id, plate_map$source)

write_csv(df,"individual_metadata.csv")

df |> dplyr::group_by(treatment) |> dplyr::summarize(ct_sex = count(sex))
