# Format_downloaded_data.R

# This script formats and prepares three files:
# - Creates "CLIN.txt"(dimension is 33-11).
# - Creates "expr_tpm.txt.gz" based on TPM data(dimension is 60591-33).
# - Creates "expr_fpkm.txt.gz" based on FPKM data(dimension is 60591-33).

# Define an input directory where all files exist
input_dir <- "~/BHK lab/kevin Project/Cholangiocarcinoma/cca_subset/"
output_dir <- "~/BHK lab/kevin Project/Cholangiocarcinoma/files/"

# Read the clinical data from the meta_data.csv file.
clinical <- read.csv(paste0(input_dir, "meta_data.csv"), stringsAsFactors = FALSE)

# Save clinical data as CLIN.txt
write.table(clinical, file = file.path(output_dir, 'CLIN.txt'), quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

# Select .txt files
file.names <- list.files(input_dir, pattern = "\\.txt$", full.names = TRUE)

# Function to process each file
process_file <- function(file_path, unique_Genename) {
  data <- read.csv(file_path, stringsAsFactors = FALSE, sep = "\t")
  data <- data[order(data$Gene.Name), ]
  data <- data[!duplicated(data$Gene.Name), ]  # Removed 35 duplicate Gene.Names.
  
  # Check if all 'Gene.Name' values are equal to 'unique_Genename'
  if (!all(data$Gene.Name == unique_Genename)) {
    stop("Gene names in ", file_path, " do not match the unique gene names.")
  }
  
  return(data)
}

# Read the first file to get the unique gene names
first_file_path <- file.names[1]  ##path to BTC_0004..txt file.
first_file_data <- read.csv(first_file_path, stringsAsFactors = FALSE, sep = "\t")
unique_Genename <- sort(unique(first_file_data$Gene.Name))

# Initialize the expr_tpm and expr_fpkm data frames
expr_tpm <- data.frame(matrix(NA, nrow = length(unique_Genename), ncol = length(clinical$id), dimnames = list(unique_Genename, clinical$id)))
expr_fpkm <- data.frame(matrix(NA, nrow = length(unique_Genename), ncol = length(clinical$id), dimnames = list(unique_Genename, clinical$id)))

# Process each file and add it to the expr_tpm and expr_fpkm data frames
for (file_path in file.names) {
  sample_data <- process_file(file_path, unique_Genename)
  # extract the sample identifier from the file name eg. "BTC_xxxx"
  sample_id <- sub("^(BTC_\\d{4}).*$", "\\1", basename(file_path))
  expr_tpm[[sample_id]] <- sample_data$TPM
  expr_fpkm[[sample_id]] <- sample_data$FPKM
}

# Save expr_tpm as "EXPR_TPM.txt.gz".
gz_tpm <- gzfile(file.path(output_dir, 'EXPR_TPM.txt.gz'), "w")
write.table(expr_tpm, file = gz_tpm, sep = "\t", row.names = TRUE, quote = FALSE)
close(gz_tpm)

# Save expr_fpkm as "EXPR_FPKM.txt.gz".
gz_fpkm <- gzfile(file.path(output_dir, 'EXPR_FPKM.txt.gz'), "w")
write.table(expr_fpkm, file = gz_fpkm, sep = "\t", row.names = TRUE, quote = FALSE)
close(gz_fpkm)


