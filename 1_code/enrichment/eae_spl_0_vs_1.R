rm(list = setdiff(ls(), "so_all"))
setwd("/singerlab/linglin/Th17_single_cell_eae_ut")

library(Seurat)
library(dplyr)
library(tidyr)
library(fgsea)
library(Matrix)
library(openxlsx)
source("1_code/utils.R")
source("1_code/enrichment/utils.R")

# configure output directories
cargs <- commandArgs(trailingOnly = TRUE)
if (length(cargs) == 0) {
  today <- "2020-05-21"
  de_date <- "2020-04-01"
} else {
  today <- cargs[1]
  de_date <- cargs[2]
}

dir_out <- paste0("2_pipeline/enrichment/EAE_SPL_0_vs_1/", today, "/")
if (!dir.exists(dir_out)) {
  dir.create(dir_out, recursive = T)
}

### load data
## DE genes
de_results <- readRDS(paste0("2_pipeline/differential_expression/SPL_EAE_0_1/", de_date, "/results.rds"))

## DB
sigs_raw_gaublomme <- read.table("0_data/gene_lists/all_Gaublomme_2016_signatures.pgf.txt", stringsAsFactors = F)
sigs_raw_gaublomme <- sigs_raw_gaublomme %>% mutate(sig_name = paste(V1, V2, sep = "-")) %>% group_by(sig_name) %>% summarise(genes = paste(V3, collapse = ","))
signature_list_1 <- strsplit(sigs_raw_gaublomme$genes, split = ",")
names(signature_list_1) <- sigs_raw_gaublomme$sig_name
# proliferating from our clustering
sigs_proliferating_raw <- readRDS("2_pipeline/differential_expression/Cluster_UT_GFPpos_inter/2020-04-16/results.rds")
tmp <- sigs_proliferating_raw[["7"]]### cluster 7 is the proliferating cluster
sigs_proliferating <- list("Proliferating" = rownames(tmp[tmp$logFC > log2(1.5) & tmp$FDR < 0.05, ]))
# interferon (Alex)
sheet_names <- getSheetNames("0_data/gene_lists/signatures/Interferon_signatures.xlsx")
sigs_raw_interferon <- lapply(sheet_names, function(x){read.xlsx("0_data/gene_lists/signatures/Interferon_signatures.xlsx", sheet = x)})
signature_list_2 <- lapply(sigs_raw_interferon, function(x){as.character(x[,1])})
names(signature_list_2) <- sheet_names
# pathogenic th17
sigs_patho17 <- list("Pathogenic_Th17" = c("Cxcl3", "Il22", "Il3", "Ccl4", "Gzmb", "Lrmp", "Ccl5", "Casp1", "Csf2", "Ccl3", "Tbx21", "lcos", "ll7r", "Stat4", "Lgals3", "Lag3"))

signature_list <- c(signature_list_1, signature_list_2, sigs_proliferating, sigs_patho17)

#### fgsea
res <- de_results$table
test_out <- test_genesets(signature_list, res, reverse = TRUE)
wb <- createWorkbook()
addWorksheet(wb, "GSEA")
writeData(wb, sheet = "GSEA", test_out[order(test_out$padj),], rowNames = F, colNames = T)
saveWorkbook(wb, file = paste0(dir_out, "GSEA_signatures.xlsx"), overwrite = T)
saveRDS(test_out, file = paste0(dir_out, "GSEA_outputs.rds"))


test_out$pathway[test_out$padj < 0.05]
