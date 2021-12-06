### Gene set enrichment analysis for each tissue vs. SPL (UT, GFP+)

rm(list = ls())
setwd("/singerlab/linglin/Th17_single_cell_eae_ut")

library(Seurat)
library(dplyr)
library(tidyr)
library(fgsea)
library(openxlsx)
source("1_code/enrichment/utils.R")

# configure output directories
cargs <- commandArgs(trailingOnly = TRUE)
if (length(cargs) == 0) {
  today <- Sys.Date()
  de_date <- "2020-03-30"
} else {
  today <- cargs[1]
  de_date <- cargs[2]
}

dir_out <- paste0("2_pipeline/enrichment/UT_GFPpos_vs_SPL/", today, "/")
if (!dir.exists(dir_out)) {
  dir.create(dir_out, recursive = T)
}

### load data
## DE results
de_results <- readRDS(paste0("2_pipeline/differential_expression/UT_GFPpos_vs_SPL/", de_date, "/results.rds"))
## GO and KEGG
kegg <- read_parse_db("0_data/other/KEGG_2019_Mouse.txt")
go_bp <- read_parse_db("0_data/other/GO_Biological_Process_2018.txt")
go_bp <- human_to_mouse(go_bp) %>% lapply(., toupper)
go_cc <- read_parse_db("0_data/other/GO_Cellular_Component_2018.txt")
go_cc <- human_to_mouse(go_cc) %>% lapply(., toupper)
go_mf <- read_parse_db("0_data/other/GO_Molecular_Function_2018.txt")
go_mf <- human_to_mouse(go_mf) %>% lapply(., toupper)


#### fgsea
tissue_vec <- c("MLN", "PP", "SI", "COL")
gsea <- list()
wb <- createWorkbook()
gsea$kegg <- lapply(tissue_vec, function(tissue) {
  cat("KEGG", tissue, "\n")
  test_out <- test_genesets(kegg, de_results[[tissue]])
  addWorksheet(wb, sheetName = tissue)
  writeData(wb, sheet = tissue, test_out[order(test_out$padj),], rowNames = F, colNames = T)
  return(test_out)
})
saveWorkbook(wb, file = paste0(dir_out, "KEGG_all.xlsx"), overwrite = T)

wb <- createWorkbook()
gsea$go_bp <- lapply(tissue_vec, function(tissue) {
  cat("GO_Biological_Process", tissue, "\n")
  test_out <- test_genesets(go_bp, de_results[[tissue]])
  addWorksheet(wb, sheetName = tissue)
  writeData(wb, sheet = tissue, test_out[order(test_out$padj),], rowNames = F, colNames = T)
  return(test_out)
})
saveWorkbook(wb, file = paste0(dir_out, "GO_Biological_Process_all.xlsx"), overwrite = T)


wb <- createWorkbook()
gsea$go_cc <- lapply(tissue_vec, function(tissue) {
  cat("GO_Cellular_Component", tissue, "\n")
  test_out <- test_genesets(go_cc, de_results[[tissue]])
  addWorksheet(wb, sheetName = tissue)
  writeData(wb, sheet = tissue, test_out[order(test_out$padj),], rowNames = F, colNames = T)
  return(test_out)
})
saveWorkbook(wb, file = paste0(dir_out, "GO_Cellular_Component_all.xlsx"), overwrite = T)

wb <- createWorkbook()
gsea$go_mf <- lapply(tissue_vec, function(tissue) {
  cat("GO_Molecular_Function", tissue, "\n")
  test_out <- test_genesets(go_mf, de_results[[tissue]])
  addWorksheet(wb, sheetName = tissue)
  writeData(wb, sheet = tissue, test_out[order(test_out$padj),], rowNames = F, colNames = T)
  return(test_out)
})
saveWorkbook(wb, file = paste0(dir_out, "GO_Molecular_Function_all.xlsx"), overwrite = T)

## save results from all DB as RDS object
saveRDS(gsea, file = paste0(dir_out, "GSEA_outputs.rds"))

