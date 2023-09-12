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
info = read_csv('/Genomics/ayroleslab2/scott/git/chromium/individual_metadata.csv')
info$name = gsub('2428__','',info$indiv)

colnames(alt_mat) <- info$name
info$group = paste0(info$plate,"_",info$treatment)

# Calculate Allele Frequency
calculate_AF <- function(alt, ref) {
  sum_alt <- sum(alt, na.rm = TRUE)
  sum_ref <- sum(ref, na.rm = TRUE)
  if ((sum_alt + sum_ref) == 0) return(NA)  # Handle cases where both counts are zero or NA
  return(sum_alt / (sum_alt + sum_ref))
}

# Function to get allele frequency for a specific site and group
get_AF_for_site_group <- function(site_index, group) {
  inds <- which(info$group == group)  # Indices of individuals in the group
  
  alt_counts <- alt_mat[site_index, inds]
  ref_counts <- ref_mat[site_index, inds]
  
  return(calculate_AF(alt = alt_counts, ref = ref_counts))
}

# Calculate allele frequencies for each combination of site and group
groups <- unique(info$group)
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
# PCA on the matrix

pca <- prcomp(AF_matrix, scale = FALSE)

library(ggfortify)
# Plot the results
png('/Genomics/ayroleslab2/scott/git/chromium/analysis/figures/pca.png')

autoplot(pca)
dev.off()

