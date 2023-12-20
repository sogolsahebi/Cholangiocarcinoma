# Format_downloaded_data.R

# Script processes expr data:
# - Generates EXPR_TPM.csv with Log2(TPM+1)
# - Creates SummarizedExperiment (SE) object

# Load libraries
library(data.table)
library(SummarizedExperiment)
library(GenomicRanges)

# Directory for input/output
directory <- "~/BHK lab/kevin Project/Cholangiocarcinoma/"

# Reading and initial processing of expr_tpm data
expr_tpm <- fread(paste0(directory, "files/EXPR_TPM.txt.gz"), sep = "\t", dec = ",", stringsAsFactors = FALSE)
rownames_expr_tpm <- expr_tpm$V1
expr_tpm <- expr_tpm[,-1]

# Converting values to numeric
expr_tpm <- apply(expr_tpm, 2, as.numeric)
rownames(expr_tpm) <- rownames_expr_tpm

range(expr_tpm) # Range is 0.0 to 196317.7

# Applying Log2 transformation
expr_tpm <- log2(expr_tpm + 1)  #Range is 0.00000  to 17.58284 

# Save processed data as CSV
write.csv(expr_tpm, paste0(directory, "files/EXPR_TPM.csv"), row.names = TRUE)

# Read clinical data
clin <- fread(paste0(directory, "files/CLIN.txt"), sep = "\t", stringsAsFactors = FALSE)
stopifnot(all(colnames(expr_tpm) == clin$id))  # Ensure column name match

# Load and process GENCODE annotation
load("~/BHK lab/Annotation/Gencode.v40.annotation (8).RData")
features_df <- features_gene
assay <- expr_tpm[rownames(expr_tpm) %in% features_df$gene_name, ]

# Replace assay gene names with gene ids
gene_ids <- unlist(lapply(rownames(assay), function(assay_row){
  vals <- rownames(features_df[features_df$gene_name == assay_row, ])
  if(length(vals) > 1){
    return(vals[1])
  } else {
    return(vals)
  }
}))

rownames(assay) <- gene_ids
assay <- assay[order(rownames(assay)), ]

# Prepare rowData for SummarizedExperiment
assay_genes <- features_df[rownames(features_df) %in% gene_ids, ]
assay_genes$gene_id_no_ver <- gsub("\\..*$", "", assay_genes$gene_id)
assay_genes <- assay_genes[!is.na(assay_genes$start), ]
row_ranges <- makeGRangesFromDataFrame(assay_genes, keep.extra.columns = TRUE)

names(row_ranges) <- row_ranges$rownames

assay_list <- list()
assay_list[["expr_tpm"]] <- assay

# Create SE object
SE_CCA <- SummarizedExperiment(assays = assay_list, colData = clin, rowRanges = row_ranges)

#Save the SE experiment ac rds file..
saveRDS(SE_expr_tpm, paste0(directory, "data/SE_CCA.rds"))

