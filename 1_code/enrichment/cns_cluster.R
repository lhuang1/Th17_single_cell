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
  today <- "2020-06-04"
  de_date <- "2020-06-01"
} else {
  today <- cargs[1]
  de_date <- cargs[2]
}

dir_out <- paste0("2_pipeline/enrichment/CNS_clusters/", today, "/")
if (!dir.exists(dir_out)) {
  dir.create(dir_out, recursive = T)
}

### load data
## DE genes
de_results <- readRDS(paste0("2_pipeline/differential_expression/Cluster_EAE_GFPall_intra/", de_date, "/CNS_results.rds"))
cluster_vec <- names(de_results)

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
gsea$go_bp <- list()
wb <- createWorkbook()
for (cluster in cluster_vec) {
  cat("\tRunning GO_Biological_Process for", cluster, "\n")
  test_out <- test_genesets(go_bp, de_results[[cluster]])
  gsea$go_bp[[cluster]] <- test_out
  addWorksheet(wb, cluster)
  writeData(wb, sheet = cluster, test_out[order(test_out$padj),], rowNames = F, colNames = T)
}
saveWorkbook(wb, file = paste0(dir_out, "GO_Biological_Process_all.xlsx"), overwrite = T)


cat("KEGG\n")
gsea$kegg <- list()
wb <- createWorkbook()
for (cluster in cluster_vec) {
  cat("\tRunning KEGG for", cluster, "\n")
  test_out <- test_genesets(kegg, de_results[[cluster]])
  gsea$kegg[[cluster]] <- test_out
  addWorksheet(wb, cluster)
  writeData(wb, sheet = cluster, test_out[order(test_out$padj),], rowNames = F, colNames = T)
}
saveWorkbook(wb, file = paste0(dir_out, "KEGG_all.xlsx"), overwrite = T)


# cat("GO_Molecular_Function\n")
# gsea$go_mf <- list()
# wb <- createWorkbook()
# for (cluster in cluster_vec) {
#   cat("\tRunning GO_Molecular_Function for", cluster, "\n")
#   test_out <- test_genesets(go_mf, de_results[[cluster]])
#   gsea$go_mf[[cluster]] <- test_out
#   addWorksheet(wb, cluster)
#   writeData(wb, sheet = cluster, test_out[order(test_out$padj),], rowNames = F, colNames = T)
# }
# saveWorkbook(wb, file = paste0(dir_out, "GO_Molecular_Function_all.xlsx"), overwrite = T)
# 
# 
# cat("GO_Cellular_Component\n")
# gsea$go_cc <- list()
# wb <- createWorkbook()
# for (cluster in cluster_vec) {
#   cat("\tRunning GO_Cellular_Component for", cluster, "\n")
#   test_out <- test_genesets(go_cc, de_results[[cluster]])
#   gsea$go_cc[[cluster]] <- test_out
#   addWorksheet(wb, cluster)
#   writeData(wb, sheet = cluster, test_out[order(test_out$padj),], rowNames = F, colNames = T)
# }
# saveWorkbook(wb, file = paste0(dir_out, "GO_Cellular_Component_all.xlsx"), overwrite = T)

saveRDS(gsea, file = paste0(dir_out, "GSEA_outputs.rds"))

