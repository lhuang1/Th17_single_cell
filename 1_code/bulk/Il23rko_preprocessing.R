#### Preprocess Il23r KO bulk RNA-seq data #####

#### configuration ####
rm(list = ls())
setwd("/singerlab/linglin/Th17_single_cell_eae_ut")

today <- Sys.Date()
dir_out <- paste0("2_pipeline/Il23rko/")
if (!dir.exists(dir_out)) dir.create(dir_out, recursive = T)

library(edgeR)
library(dplyr)
library(openxlsx)
library(tibble)
library(ggplot2)


#### load data
tpm1 <- read.table(gzfile(paste0("0_data/counts/bulk/Il23rko/P1.dge.txt.gz")), header = T, row.names = 1)
tpm2 <- read.table(gzfile(paste0("0_data/counts/bulk/Il23rko/P2.dge.txt.gz")), header = T, row.names = 1)


# marker_vec <- c("Cxcr6", "Slam", "Tdt")
# replica_vec <- 1:3 ## only use mouse 2,3,4 (rename to 1,2,3 to avoid confusion); mouse 5 is UT; exclude mouse 1 because Cxcr6 has low number of reads and alignment/unique rate.
# tpm <- list()
# for (m in marker_vec) {
#   for (r in replica_vec) {
#     sample_name <- paste("SPL", m, r, sep = "_")
#     cat("Reading", sample_name, "...\n")
#     x = read.table(gzfile(paste0("0_data/counts/bulk/counts/", paste("SPL", m, r+1, sep = "_"), ".dge.txt.gz")),
#                    header = T, row.names = 1)  
#     tpm[[sample_name]] <- rowSums(x)
#   }
# }

# tpm_df <- as.data.frame(tpm)
# saveRDS(tpm_df, paste0(dir_out, "TPM.rds"))
# meta_data <- Reduce(rbind, strsplit(colnames(tpm_df), split = "_")) %>% data.frame()
# rownames(meta_data) <- colnames(tpm_df)
# colnames(meta_data) <- c("tissue", "marker", "replica")
# meta_data$marker <- factor(meta_data$marker, levels = c("Tdt", "Cxcr6", "Slam"))
# 
# ##### Use Tdt as reference #####
# # Create dge object
# dgeFull <- DGEList(tpm_df)
# # Estimate the normalization factors
# dgeFull <- calcNormFactors(dgeFull, method="TMM")
# # create design matrix
# # design <- model.matrix(~ marker + replica, data = meta_data)
# design <- model.matrix(~ marker, data = meta_data) # ignore pairing
# # Estimate dispersion
# dgeFull <- estimateDisp(dgeFull, design = design)
# # Perform QLF tests
# fit <- glmFit(dgeFull, design = design)
# lrt_Cxcr6 <- glmLRT(fit, coef = 2)
# lrt_Slam <- glmLRT(fit, coef = 3)
# res_Cxcr6 <- topTags(lrt_Cxcr6, n = Inf, adjust.method = "BH", sort.by = "none", p.value = 1)
# res_Slam <- topTags(lrt_Slam, n = Inf, adjust.method = "BH", sort.by = "none", p.value = 1)
# # save DE results
# saveRDS(res_Cxcr6, file = paste0(dir_out, "edger_cxcr6_vs_tdt.rds"))
# saveRDS(res_Slam, file = paste0(dir_out, "edger_slam_vs_tdt.rds"))
# 
# ## diagnostic plot
# pdf(paste0(dir_out, "bulk_de_diagnostic_vs_tdt.pdf"))
# df_cxcr6_vs_tdt <- data.frame(res_Cxcr6$table)
# ggplot(df_cxcr6_vs_tdt) + 
#   geom_point(aes(x = logCPM, y = logFC, color = FDR < 0.05)) +
#   scale_color_manual(values = c("black", "red")) +
#   ggtitle("Cxcr6")
# ggplot(df_cxcr6_vs_tdt) + 
#   geom_point(aes(x = logFC, y = -log10(PValue), color = FDR < 0.05)) +
#   scale_color_manual(values = c("black", "red")) +
#   ggtitle("Cxcr6")
# df_slamf6_vs_tdt <- data.frame(res_Slam$table)
# ggplot(df_slamf6_vs_tdt) + 
#   geom_point(aes(x = logCPM, y = logFC, color = FDR < 0.05)) +
#   scale_color_manual(values = c("black", "red")) +
#   ggtitle("Slamf6")
# ggplot(df_slamf6_vs_tdt) + 
#   geom_point(aes(x = logFC, y = -log10(PValue), color = FDR < 0.05)) +
#   scale_color_manual(values = c("black", "red")) +
#   ggtitle("Slamf6")
# dev.off()
# ## obtain signatures for each cluster
# sigs_vs_tdt <- list()
# sigs_vs_tdt[["Slamf6-plus"]] <- df_slamf6_vs_tdt %>% 
#   rownames_to_column() %>% 
#   filter(FDR < 0.05 & logFC > log2(1.5)) %>% 
#   arrange(FDR) %>%  ## sort by FDR increasing order
#   select(rowname) %>% 
#   unlist() %>% 
#   as.character()
# sigs_vs_tdt[["Slamf6-minus"]] <- df_slamf6_vs_tdt %>% 
#   rownames_to_column() %>% 
#   filter(FDR < 0.05 & logFC < -log2(1.5)) %>% 
#   arrange(FDR) %>% 
#   select(rowname) %>% 
#   unlist() %>% 
#   as.character()
# sigs_vs_tdt[["Cxcr6-plus"]] <- df_cxcr6_vs_tdt %>% 
#   rownames_to_column() %>% 
#   filter(FDR < 0.05 & logFC > log2(1.5)) %>% 
#   arrange(FDR) %>% 
#   select(rowname) %>% 
#   unlist() %>% 
#   as.character()
# sigs_vs_tdt[["Cxcr6-minus"]] <- df_cxcr6_vs_tdt %>% 
#   rownames_to_column() %>% 
#   filter(FDR < 0.05 & logFC < -log2(1.5)) %>% 
#   arrange(FDR) %>% 
#   select(rowname) %>% 
#   unlist() %>% 
#   as.character()
# 
# wb <- createWorkbook()
# sapply(names(sigs_vs_tdt), function(i){
#   addWorksheet(wb, sheetName = i)
#   writeData(wb, sheet = i, x = sigs_vs_tdt[[i]])
#   length(sigs_vs_tdt[[i]])
# })
# saveWorkbook(wb, file = paste0(dir_out, "signature_vs_tdt_all.xlsx"), overwrite = T)
# saveRDS(sigs_vs_tdt, file = paste0(dir_out, "signature_vs_tdt_all.rds"))
# 
# ## filter out microRNA and small non-coding RNAs
# wb <- createWorkbook()
# signature_vs_tdt_filtered <- sapply(names(sigs_vs_tdt), function(i){
#   tmp = sigs_vs_tdt[[i]]
#   tmp = tmp[!(grepl("Mir", tmp) | grepl("Snord", tmp))]
#   addWorksheet(wb, sheetName = i)
#   writeData(wb, sheet = i, x = tmp)
#   cat(i, ": ", length(tmp), "\n")
#   tmp
# })
# saveWorkbook(wb, file = paste0(dir_out, "signature_vs_tdt_filtered.xlsx"), overwrite = T)
# saveRDS(signature_vs_tdt_filtered, file = paste0(dir_out, "signature_vs_tdt_filtered.rds"))
# 
# 
# 
# ##### Direct comparison between cxcr6 and slamf6 #####
# is_tdt <- grepl("Tdt", colnames(tpm_df))
# # Create dge object
# dge_sub <- DGEList(tpm_df[!is_tdt])
# # Estimate the normalization factors
# dge_sub <- calcNormFactors(dge_sub, method="TMM")
# # create design matrix
# # design <- model.matrix(~ marker + replica, data = droplevels(meta_data[!is_tdt,]))
# design <- model.matrix(~ marker, data = droplevels(meta_data[!is_tdt,])) # ignore pairing
# # Estimate dispersion
# dge_sub <- estimateDisp(dge_sub, design = design)
# # Perform QLF tests
# fit <- glmFit(dge_sub, design = design)
# lrt_Slam_vs_Cxcr6 <- glmLRT(fit, coef = 2)
# res_Slam_vs_Cxcr6 <- topTags(lrt_Slam_vs_Cxcr6, n = Inf, adjust.method = "BH", sort.by = "none", p.value = 1)
# # save DE results
# saveRDS(res_Slam_vs_Cxcr6, file = paste0(dir_out, "edger_slam_vs_cxcr6.rds"))
# 
# ## diagnostic plot
# pdf(paste0(dir_out, "bulk_de_diagnostic_slamf6_vs_cxcr6.pdf"))
# df_slamf6_vs_cxcr6 <- data.frame(res_Slam_vs_Cxcr6$table)
# ggplot(df_slamf6_vs_cxcr6) + 
#   geom_point(aes(x = logCPM, y = logFC, color = FDR < 0.05)) +
#   scale_color_manual(values = c("black", "red")) +
#   ggtitle("Slamf6 vs. Cxcr6")
# ggplot(df_cxcr6_vs_tdt) + 
#   geom_point(aes(x = logFC, y = -log10(PValue), color = FDR < 0.05)) +
#   scale_color_manual(values = c("black", "red")) +
#   ggtitle("Slamf6 vs. Cxcr6")
# dev.off()
# ## obtain signatures for each cluster
# sigs_slamf6_vs_cxcr6 <- list()
# sigs_slamf6_vs_cxcr6[["Slamf6-plus"]] <- df_slamf6_vs_cxcr6 %>% 
#   rownames_to_column() %>% 
#   filter(FDR < 0.05 & logFC > log2(1.5)) %>% 
#   arrange(FDR) %>%  ## sort by FDR increasing order
#   select(rowname) %>% 
#   unlist() %>% 
#   as.character()
# sigs_slamf6_vs_cxcr6[["Cxcr6-plus"]] <- df_slamf6_vs_cxcr6 %>% 
#   rownames_to_column() %>% 
#   filter(FDR < 0.05 & logFC < -log2(1.5)) %>% 
#   arrange(FDR) %>% 
#   select(rowname) %>% 
#   unlist() %>% 
#   as.character()
# 
# wb <- createWorkbook()
# sapply(names(sigs_slamf6_vs_cxcr6), function(i){
#   addWorksheet(wb, sheetName = i)
#   writeData(wb, sheet = i, x = sigs_slamf6_vs_cxcr6[[i]])
#   cat(i, ": ", length(sigs_slamf6_vs_cxcr6[[i]]), "\n")
# })
# saveWorkbook(wb, file = paste0(dir_out, "signature_slamf6_vs_cxcr6_all.xlsx"), overwrite = T)
# saveRDS(sigs_slamf6_vs_cxcr6, file = paste0(dir_out, "signature_slamf6_vs_cxcr6_all.rds"))
# 
# ## filter out microRNA and small non-coding RNAs
# wb <- createWorkbook()
# sigs_slamf6_vs_cxcr6_filtered <- sapply(names(sigs_slamf6_vs_cxcr6), function(i){
#   tmp = sigs_slamf6_vs_cxcr6[[i]]
#   tmp = tmp[!(grepl("Mir", tmp) | grepl("Snord", tmp))]
#   addWorksheet(wb, sheetName = i)
#   writeData(wb, sheet = i, x = tmp)
#   cat(i, ": ", length(tmp), "\n")
#   tmp
# })
# saveWorkbook(wb, file = paste0(dir_out, "signature_slamf6_vs_cxcr6_filtered.xlsx"), overwrite = T)
# saveRDS(sigs_slamf6_vs_cxcr6_filtered, file = paste0(dir_out, "signature_slamf6_vs_cxcr6_filtered.rds"))
