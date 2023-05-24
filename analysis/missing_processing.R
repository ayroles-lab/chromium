# Load the IDs list from a file
id_file <- "/Genomics/ayroleslab2/scott/git/chromium/sequence_processing/list_of_samples.txt"  # Replace with the actual path to your IDs file
id_list <- readLines(id_file)

# Specify the folder to search for files
folder_path <- "path/to/folder"  # Replace with the actual path to the folder containing the files

# Get a list of all files in the folder
files <- list.files("/Genomics/ayroleslab2/scott/git/chromium/data/raw/", full.names = TRUE)

# Filter the files to retain only the .vcf.gz.tbi files
vcf_files <- files[grepl("\\.vcf\\.gz\\.tbi$", files)]

# Extract the IDs from the file names
file_ids <- gsub("\\.recal\\.vcf\\.gz\\.tbi$", "", basename(vcf_files))

# Find the IDs that don't have matching files
missing_ids <- setdiff(id_list, file_ids)

# Find the locations of missing IDs in the list
missing_ids_location <- which(id_list %in% missing_ids)

missing_df <- data.frame(ID = missing_ids, Location = missing_ids_location)

# Print the data frame
print(missing_df)

paste(missing_df$Location, collapse = ",")
