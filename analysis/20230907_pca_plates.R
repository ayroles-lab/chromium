library(data.table)
library(tidyverse)
library(parallel)

# Define paths and list files
path <- "/Genomics/ayroleslab2/scott/git/chromium/data/snakemake_tmp/count_info/archive/"
alt_count_files <- list.files(path=path, pattern="alternate_counts*", full.names=TRUE)
ref_count_files <- list.files(path=path, pattern="reference_counts_*", full.names=TRUE)

# Read and merge alternate and reference count data
alt_count_data <- rbindlist(lapply(alt_count_files, fread), fill=TRUE)
ref_count_data <- rbindlist(lapply(ref_count_files, fread), fill=TRUE)


# Convert them to matrices for faster computation
alt_mat <- as.matrix(alt_count_data[, -1, with=FALSE])
ref_mat <- as.matrix(ref_count_data[, -1, with=FALSE])

row.names(alt_mat) <- alt_count_data$site
row.names(ref_mat) <- ref_count_data$site


# Read the info file
info = read_csv('/Genomics/ayroleslab2/scott/git/chromium/metadata/longevity_dna_individual_metadata.csv')
info$name = gsub('2428__','',info$indiv)

colnames(alt_mat) <- info$name
colnames(ref_mat) <- info$name


# -------------------------------

# Drop plates
filtered_info <- info %>% filter(!(plate %in% c("Chrom_31_7","Chrom_31_8","Chrom_31_9")))

ref_mat <- ref_mat[, filtered_info$name]
alt_mat <- alt_mat[, filtered_info$name]

info <- filtered_info


colnames(alt_mat) <- info$name
info$group = paste0(info$cage,"_",info$treatment)

# Calculate Allele Frequency
calculate_AF <- function(alt, ref) {
  sum_alt <- sum(alt, na.rm = TRUE)
  sum_ref <- sum(ref, na.rm = TRUE)
  if ((sum_alt + sum_ref) == 0) return(NA)  # Handle cases where both counts are zero or NA
  return(sum_alt / (sum_alt + sum_ref))
}

# Function to get allele frequency for a specific site and group
get_AF_for_site_group <- function(site_index, group) {
  inds <- which(info$plate == group)  # Indices of individuals in the group
  
  alt_counts <- alt_mat[site_index, inds]
  ref_counts <- ref_mat[site_index, inds]

  # Randomly draw allele for each individual
  p = alt_counts / (alt_counts + ref_counts)
  drawn_alleles = rbinom(length(p), 1, p)
  
  # Use the drawn alleles to compute allele frequency
  return(calculate_AF(alt = drawn_alleles, ref = 1-drawn_alleles))
}

# Calculate allele frequencies for each combination of site and group
groups <- unique(info$plate)
sites <- row.names(ref_mat)

# Using expand.grid to get every combination of site and group
combinations <- expand.grid(site = sites, group = groups)


# Get number of cores
num_cores <- 48

# Function to compute allele frequencies for a chunk of combinations
get_AF_for_chunk <- function(chunk) {
  mapply(get_AF_for_site_group, 
         site_index = chunk$site, 
         group = chunk$group)
}

# Split combinations into chunks for parallel processing
chunks <- split(combinations, ceiling(seq_len(nrow(combinations)) / (nrow(combinations) / num_cores)))

# Use mclapply to process each chunk in parallel
results_list <- mclapply(chunks, get_AF_for_chunk, mc.cores = num_cores)

combinations$allele_frequency <- unlist(results_list)

# Convert to a data.table
result <- data.table(combinations)
wide_result <- result %>%
  pivot_wider(names_from = group, values_from = allele_frequency)

# Convert to a matrix
AF_matrix <- as.matrix(wide_result[-1])
rownames(AF_matrix) <- wide_result$site

# Drop rows with NAs
AF_matrix <- AF_matrix[complete.cases(AF_matrix), ]

# Run PCA on the transposed AF_matrix
pca_result <- prcomp(t(AF_matrix))

# Convert PCA results to a data frame for ggplot
pca_df <- as.data.frame(pca_result$x)

# Add groups to the data frame (this time, directly extracted from the column names)
pca_df$group <- factor(colnames(AF_matrix))
df <- data.frame(x = colnames(AF_matrix))
df <- df %>% separate(x,sep="_", into = c("cage", "treatment"))
pca_df <- cbind(pca_df, df)
pve <- round(pca_result$sdev^2 / sum(pca_result$sdev^2) * 100, 1)

# Plot the first two principal components
ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  labs(title = "PCA of Allele Frequencies -- Grouped by Plate", x = paste0("PC1 (", pve[1], "%)"), y = paste0("PC2 (", pve[2], "%)")) +
  theme(legend.position = "bottom") +
  geom_label(aes(label = group), nudge_x = 0.1, nudge_y = 0.1, size = 3)

ggsave("/Genomics/ayroleslab2/scott/git/chromium/analysis/figures/pca_af_plates_689_removed.jpeg", width = 10, height = 10, units = "in")


summary(pca_result)
