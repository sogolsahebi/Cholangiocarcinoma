# Format_downloaded_data.R.

# This script formats and Prepare three files.
# - Creates "CLIN.txt".
# - Creates "EXPR_TPM.txt.gz"  based on TPM data.
# - Creates "EXPR_FPKM.txt.gz" based on FPMK data.



#Define a input_dir where all files exists.
input_dir <- "~/BHK lab/kevin Project/Cholangiocarcinoma/cca_subset/"

#1. CLIN.txt
#Read the meta_data,csv as set it as clinical data.
clinical <- read.csv(paste0(input_dir, "meta_data.csv"), stringsAsFactors = FALSE) #dimension  33-11

#TODO: save the ClIN.txt file.

#2. Preparing for EXPR_TPM file. 

#TODO: reading all --33?-- files. 

#1 BTC_0004 sample with dimention 60626- 9
BTC_0004 <- read.csv(paste0(input_dir, "BTC_0004_Lv_M_526_stringtie_abundance.txt"), stringsAsFactors = FALSE, sep = "\t")

#remove the duplication Gene.Name, keeping the first occurrence

sum(duplicated(BTC_0004$Gene.Name)) #35 duplicate Gene.Names

BTC_0004 <- BTC_0004[!duplicated(BTC_0004$Gene.Name), ] #now dimension is 60591-9

#Store The Gene.Names and set it as rownames of the data.table for EXPM_TPM.
unique_Genename <- BTC_0004$Gene.Name

#Create a data.table with Row names of Gene.Name and columnnames of out clinical id.
# Create the empty data frame with specified row and column names #60591- 33
EXPR_TPM <- data.frame(matrix(NA, nrow = length(BTC_0004$Gene.Name), ncol = length(clinical$id), 
                              dimnames = list(BTC_0004$Gene.Name, clinical$id))) 

#Set BTC_0004 clomuns of EXPR_TPM.
EXPR_TPM$BTC_0004 <- BTC_0004$TPM 

#Adding the Remining 32 files tpm data to EXPR_TPM file.
#double check if we  have same Gene.names as BTC_0004 Gene.names("unique_Genename").
BTC_0009 <- read.csv(paste0(input_dir, "BTC_0009_Lv_P_526_stringtie_abundance"), stringsAsFactors = FALSE, sep = "\t")


BTC_0010_Lv_P_526_stringtie_abundance
BTC_0011_Lv_P_526_stringtie_abundance
BTC_0012_Lv_P_526_stringtie_abundance
BTC_0013_Lv_P_526_stringtie_abundance
BTC_0014_Lv_M_526_stringtie_abundance
BTC_0015_Lv_P_526_stringtie_abundance
BTC_0016_Lv_M_526_stringtie_abundance
BTC_0017_Lv_P_526_stringtie_abundance
BTC_0018_Lv_P_526_stringtie_abundance
BTC_0019_Lv_P_526_stringtie_abundance
BTC_0020_Ab_M_526_stringtie_abundance
BTC_0022_Lv_M_526_stringtie_abundance
BTC_0023_Lv_M_526_stringtie_abundance
BTC_0024_Lv_M_526_stringtie_abundance
BTC_0025_Sn_M_526_stringtie_abundance
BTC_0026_Lv_M_526_stringtie_abundance
BTC_0028_Lv_M_526_stringtie_abundance
BTC_0029_Lv_M_526_stringtie_abundance
BTC_0030_Lv_M_526_stringtie_abundance
BTC_0032_Lv_M_526_stringtie_abundance
BTC_0033_Lv_P_526_stringtie_abundance
BTC_0034_Lv_P_526_stringtie_abundance
BTC_0035_Lv_M_526_stringtie_abundance
BTC_0036_Lv_P_526_stringtie_abundance
BTC_0037_Lv_P_526_stringtie_abundance
BTC_0038_Lv_P_526_stringtie_abundance
BTC_0039_Lv_M_526_stringtie_abundance
BTC_0040_Lv_P_526_stringtie_abundance
BTC_0042_Lv_P_526_stringtie_abundance
BTC_0045_Lv_P_526_stringtie_abundance
BTC_8002_Lv_M_526_stringtie_abundance

