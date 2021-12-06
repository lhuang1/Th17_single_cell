rm(list = ls())
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
  today <- Sys.Date()
  de_date <- "2020-05-11"
} else {
  today <- cargs[1]
  de_date <- cargs[2]
}

dir_out <- paste0("2_pipeline/enrichment/EAE_SI_5_cluster/", today, "/")
if (!dir.exists(dir_out)) {
  dir.create(dir_out, recursive = T)
}

### load data
## DE genes
de_results <- readRDS(paste0("2_pipeline/differential_expression/Cluster_EAE_GFPall_intra/", de_date, "/SI_5_results.rds"))

## DB
go_bp <- read_parse_db("0_data/other/GO_Biological_Process_2018.txt")
go_bp <- human_to_mouse(go_bp) %>% lapply(., toupper)
go_mf <- read_parse_db("0_data/other/GO_Molecular_Function_2018.txt")
go_mf <- human_to_mouse(go_mf) %>% lapply(., toupper)
go_cc <- read_parse_db("0_data/other/GO_Cellular_Component_2018.txt")
go_cc <- human_to_mouse(go_cc) %>% lapply(., toupper)
kegg <- read_parse_db("0_data/other/KEGG_2019_Mouse.txt")

#### fgsea
gsea <- list()

cat("GO_Biological_Process\n")
test_out <- test_genesets(go_bp, de_results)
gsea$go_bp <- test_out
wb <- createWorkbook()
addWorksheet(wb, "GSEA")
writeData(wb, sheet = "GSEA", test_out[order(test_out$padj),], rowNames = F, colNames = T)
saveWorkbook(wb, file = paste0(dir_out, "GO_Biological_Process_all.xlsx"), overwrite = T)


cat("KEGG\n")
test_out <- test_genesets(kegg, de_results)
gsea$kegg <- test_out
wb <- createWorkbook()
addWorksheet(wb, "GSEA")
writeData(wb, sheet = "GSEA", test_out[order(test_out$padj),], rowNames = F, colNames = T)
saveWorkbook(wb, file = paste0(dir_out, "KEGG_all.xlsx"), overwrite = T)


cat("GO_Molecular_Function\n")
test_out <- test_genesets(go_mf, de_results)
gsea$go_mf <- test_out
wb <- createWorkbook()
addWorksheet(wb, "GSEA")
writeData(wb, sheet = "GSEA", test_out[order(test_out$padj),], rowNames = F, colNames = T)
saveWorkbook(wb, file = paste0(dir_out, "GO_Molecular_Function_all.xlsx"), overwrite = T)


cat("GO_Cellular_Component\n")
test_out <- test_genesets(go_cc, de_results)
gsea$go_cc <- test_out
wb <- createWorkbook()
addWorksheet(wb, "GSEA")
writeData(wb, sheet = "GSEA", test_out[order(test_out$padj),], rowNames = F, colNames = T)
saveWorkbook(wb, file = paste0(dir_out, "GO_Cellular_Component_all.xlsx"), overwrite = T)

## save results from all DB as RDS object
saveRDS(gsea, file = paste0(dir_out, "GSEA_outputs.rds"))

