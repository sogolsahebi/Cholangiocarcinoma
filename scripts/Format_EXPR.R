# Script to Process Expression Data
# - Generates EXPR_TPM.csv with Log2(TPM+1)
# - Creates SummarizedExperiment (SE) Object - annottation : used GENCODE 40 and 38040 genes was found.

# Load Required Libraries
library(data.table) 
library(SummarizedExperiment)
library(GenomicRanges)

# Set Directory for Input/Output
directory <- "~/BHK lab/kevin Project/Cholangiocarcinoma/"

# Read and Process Expression TPM Data
expr_tpm_file <- file.path(directory, "files/EXPR_TPM.txt.gz")
expr_tpm <- fread(expr_tpm_file, sep = "\t", dec = ",", stringsAsFactors = FALSE)
rownames_expr_tpm <- expr_tpm$V1
expr_tpm <- expr_tpm[,-1]

# Convert Values to Numeric and Assign Row Names
expr_tpm <- apply(expr_tpm, 2, as.numeric)
rownames(expr_tpm) <- rownames_expr_tpm

# Log2 Transformation
expr_tpm <- log2(expr_tpm + 1)

# Save Processed Data as CSV
expr_tpm_csv_file <- paste0(directory, "files/EXPR_TPM.csv")
write.csv(expr_tpm, expr_tpm_csv_file, row.names = TRUE)

# Read Clinical Data and Ensure Column Name Match
clin_file <- paste0(directory, "files/CLIN.txt")
clin <- fread(clin_file, sep = "\t", stringsAsFactors = FALSE)
stopifnot(all(colnames(expr_tpm) == clin$id))

# Load and Process GENCODE Annotation
gencode_file <- "~/BHK lab/Annotation/Gencode.v40.annotation.RData"
load(gencode_file)
features_df <- features_gene

# Remove Duplicated Gene Names
features_df <- features_df[!duplicated(features_df$gene_name), ]

# Filter and Order Assay Data
assay <- expr_tpm[rownames(expr_tpm) %in% features_df$gene_name, ]
assay <- assay[order(rownames(assay)), ]

# Prepare Row Data for SummarizedExperiment
assay_genes <- features_df[features_df$gene_name %in% rownames(assay), ]
assay_genes$gene_id_no_ver <- gsub("\\..*$", "", assay_genes$gene_id)
assay_genes <- assay_genes[!is.na(assay_genes$start), ]
rownames(assay_genes) <- assay_genes$gene_name
assay_genes <- assay_genes[order(rownames(assay_genes)), ]

# Create GRanges Object for Row Ranges
row_ranges <- makeGRangesFromDataFrame(assay_genes, keep.extra.columns = TRUE)
names(row_ranges) <- assay_genes$gene_name

# Create Assay List
assay_list <- list(expr_tpm = assay)

# Create SummarizedExperiment Object
SE_CCA <- SummarizedExperiment(assays = assay_list, colData = clin, rowRanges = row_ranges)

# Save the SummarizedExperiment Object
se_cca_file <- paste0(directory, "data/SE_CCA.rds")
saveRDS(SE_CCA, se_cca_file)



