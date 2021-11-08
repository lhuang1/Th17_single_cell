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
  today <- "2020-03-12"
  clustering_date <- "2020-03-03"
} else {
  today <- cargs[1]
  clustering_date <- cargs[2]
}

dir_out <- paste0("2_pipeline/enrichment/UT_GFPpos_vs_SPL/", today, "/")
if (!dir.exists(dir_out)) {
  dir.create(dir_out, recursive = T)
}

### load data
so_all <- readRDS("2_pipeline/preprocessing/integration/2020-01-29/so_integrated_SCTransformed.rds")
keep_cells_clean <- read.table(paste0("2_pipeline/clustering/UT_GFPpos_inter/", clustering_date, "/FILES/meta_data.txt"))
so <- so_all[,rownames(keep_cells_clean)]
so@meta.data <- droplevels(so@meta.data)
rm(so_all)

## DE genes
de_results <- readRDS("2_pipeline/differential_expression/UT_GFPpos_vs_SPL/2020-03-03/results.rds")

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


saveRDS(gsea, file = paste0(dir_out, "GSEA_outputs.rds"))


### combine GO_BP lists for different tissues
results <- readRDS(paste0(dir_out, "GSEA_outputs.rds"))
# df <- sapply(gsea$go_bp, function(x){
#   x$padj
# }) %>% data.frame()
# 
# colnames(df) <- tissue_vec
# df <- cbind(pathway = gsea$go_bp[[1]]$pathway, df)
# 
# 
# wb <- createWorkbook()
# sigStyle <- createStyle(bgFill = "#FFC7CE")
# 
# addWorksheet(wb, sheetName = "sig_terms")
# n_tissue_sig <- rowSums(df[,tissue_vec] < 0.05)
# df_sig <- df[n_tissue_sig > 0 & !is.na(n_tissue_sig),]
# writeData(wb, sheet = "sig_terms", x = df_sig, colNames = T)
# conditionalFormatting(wb, "sig_terms", cols = 1:ncol(df_sig), rows = 1:(nrow(df_sig) + 1), rule="<0.05", style = sigStyle)
# 
# addWorksheet(wb, sheetName = "all_terms")
# writeData(wb, sheet = "all_terms", x = df, colNames = T)
# conditionalFormatting(wb, "all_terms", cols = 1:ncol(df), rows = 1:(nrow(df) + 1), rule="<0.05", style = sigStyle)
# 
# saveWorkbook(wb, file = paste0(dir_out, "GO_Biological_Process_all_combined.xlsx"), overwrite = T)
# 
# library(DESeq2)

### combine GO_BP lists for different tissues
today <- "2020-12-24"
wb <- createWorkbook()
sigStyle <- createStyle(bgFill = "#FFC7CE")

## by p-value
df <- sapply(results$go_bp, function(x){
  ifelse(x$NES > 0, x$pval, 1)
}) %>% data.frame()
colnames(df) <- tissue_vec
df <- cbind(pathway = results$go_bp[[1]]$pathway, df)
n_tissue_sig <- rowSums(df[,tissue_vec] < 0.05)
df_sig <- df[n_tissue_sig > 0 & !is.na(n_tissue_sig),]
sig_value <- (df_sig[,tissue_vec] < 0.05) %*% c(1.3, 1.2, 1.1, 1.0)
df_sig <- df_sig[order(sig_value, decreasing = T),]
sheet <- "pval_sig_terms"
addWorksheet(wb, sheetName = sheet)
writeData(wb, sheet = sheet, x = df_sig, colNames = T)
conditionalFormatting(wb, sheet, cols = 1:ncol(df_sig), rows = 1:(nrow(df_sig) + 1), rule="<0.05", style = sigStyle)
sheet <- "pval_all_terms"
addWorksheet(wb, sheetName = sheet)
writeData(wb, sheet = sheet, x = df, colNames = T)
conditionalFormatting(wb, sheet, cols = 1:ncol(df), rows = 1:(nrow(df) + 1), rule="<0.05", style = sigStyle)

## by p-value
df <- sapply(results$go_bp, function(x){
  ifelse(x$NES > 0, x$padj, 1)
}) %>% data.frame()
colnames(df) <- tissue_vec
df <- cbind(pathway = results$go_bp[[1]]$pathway, df)
n_tissue_sig <- rowSums(df[,tissue_vec] < 0.05)
df_sig <- df[n_tissue_sig > 0 & !is.na(n_tissue_sig),]
sig_value <- (df_sig[,tissue_vec] < 0.05) %*% c(1.3, 1.2, 1.1, 1.0)
df_sig <- df_sig[order(sig_value, decreasing = T),]
sheet <- "fdr_sig_terms"
addWorksheet(wb, sheetName = sheet)
writeData(wb, sheet = sheet, x = df_sig, colNames = T)
conditionalFormatting(wb, sheet, cols = 1:ncol(df_sig), rows = 1:(nrow(df_sig) + 1), rule="<0.05", style = sigStyle)
sheet <- "fdr_all_terms"
addWorksheet(wb, sheetName = sheet)
writeData(wb, sheet = sheet, x = df, colNames = T)
conditionalFormatting(wb, sheet, cols = 1:ncol(df), rows = 1:(nrow(df) + 1), rule="<0.05", style = sigStyle)

saveWorkbook(wb, file = paste0(dir_out, "GSEA_GO_Biological_Process_all_combined_", today, ".xlsx"), overwrite = T)

