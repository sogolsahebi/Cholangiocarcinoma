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

# Define the unique part of the sample file names.
sample_files <- c(
  "BTC_0004_Lv_M_", "BTC_0009_Lv_P_", "BTC_0010_Lv_P_", "BTC_0011_Lv_P_", 
  "BTC_0012_Lv_P_", "BTC_0013_Lv_P_", "BTC_0014_Lv_M_", "BTC_0015_Lv_P_", 
  "BTC_0016_Lv_M_", "BTC_0017_Lv_P_", "BTC_0018_Lv_P_", "BTC_0019_Lv_P_", 
  "BTC_0020_Ab_M_", "BTC_0022_Lv_M_", "BTC_0023_Lv_M_", "BTC_0024_Lv_M_", 
  "BTC_0025_Sn_M_", "BTC_0026_Lv_M_", "BTC_0028_Lv_M_", "BTC_0029_Lv_M_", 
  "BTC_0030_Lv_M_", "BTC_0032_Lv_M_", "BTC_0033_Lv_P_", "BTC_0034_Lv_P_", 
  "BTC_0035_Lv_M_", "BTC_0036_Lv_P_", "BTC_0037_Lv_P_", "BTC_0038_Lv_P_", 
  "BTC_0039_Lv_M_", "BTC_0040_Lv_P_", "BTC_0042_Lv_P_", "BTC_0045_Lv_P_", 
  "BTC_8002_Lv_M_"
)

# Function to process each file
process_file <- function(unique_file_name, input_dir, unique_Genename) {
  file_name <- paste0(unique_file_name, "526_stringtie_abundance.txt")
  file_path <- paste0(input_dir, file_name)
  data <- read.csv(file_path, stringsAsFactors = FALSE, sep = "\t")
  data <- data[order(data$Gene.Name), ]
  data <- data[!duplicated(data$Gene.Name), ] # Removed 35 duplicate Gene.Names.
  
  # Check if all 'Gene.Name' values are equal to 'unique_Genename'
  if (!all(data$Gene.Name == unique_Genename)) {
    stop("Gene names in ", file_name, " do not match the unique gene names.")
  }
  
  return(data)
}

# Read the first file to get the unique gene names
first_file_path <- paste0(input_dir, sample_files[1], "526_stringtie_abundance.txt") #path to BTC_0004..txt file.
first_file_data <- read.csv(first_file_path, stringsAsFactors = FALSE, sep = "\t")
unique_Genename <- sort(unique(first_file_data$Gene.Name))

# Initialize the expr_tpm and expr_fpkm data frames
expr_tpm <- data.frame(matrix(NA, nrow = length(unique_Genename), ncol = length(clinical$id),
                              dimnames = list(unique_Genename, clinical$id)))

expr_fpkm <- data.frame(matrix(NA, nrow = length(unique_Genename), ncol = length(clinical$id),
                               dimnames = list(unique_Genename, clinical$id)))

# Process each file and add it to the expr_tpm and expr_fpkm data frames
for (file_name in sample_files) {
  sample_data <- process_file(file_name, input_dir, unique_Genename)
  # extract the sample identifier from the file name eg. "BTC_xxxx"
  sample_id <- sub("^(BTC_\\d{4}).*$", "\\1", file_name)
  expr_tpm[[sample_id]] <- sample_data$TPM
  expr_fpkm[[sample_id]] <- sample_data$FPKM
}

# Save expr_tpm as "expr_tpm.txt.gz".
gz_tpm <- gzfile(file.path(output_dir, 'expr_tpm.txt.gz'), "w")
write.table(expr_tpm, file = gz_tpm, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
close(gz_tpm)

#Save as csv too.
#write.csv(expr_tpm, paste0(output_dir, "EXPR_TPM.csv"), row.names = TRUE)

# Save expr_fpkm as "expr_fpkm.txt.gz".
gz_fpkm <- gzfile(file.path(output_dir, 'expr_fpkm.txt.gz'), "w")
write.table(expr_fpkm, file = gz_fpkm, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
close(gz_fpkm)

#Save as csv too.
#write.csv(expr_tpm, paste0(output_dir, "EXPR_FPKM.csv"), row.names = TRUE)
