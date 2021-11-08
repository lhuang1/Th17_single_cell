####### Analyze human (Schafflick) data ########

## configuration
rm(list = ls())
setwd("/singerlab/linglin/Th17_single_cell_eae_ut/")

library(Seurat)
library(ggplot2)
library(gridExtra)
library(ggsignif)
library(ggridges)
library(dplyr)
library(tibble)
library(openxlsx)
source("1_code/utils.R")

integrate_date <- "2020-06-01"
clustering_date <- "2020-06-11"
today <- "2020-07-01"
dir_out <- paste0("2_pipeline/human_shafflick/", today, "/")
if (!dir.exists(dir_out)) {
  dir.create(dir_out, recursive = T)
}

##### load data ######
pbmc <- readRDS(file = paste0("2_pipeline/human_shafflick/", integrate_date, "/so_combined_PBMCs.rds"))
meta <- readRDS(file = paste0("2_pipeline/human_shafflick/", clustering_date, "/meta_CD4_PBMCs.rds"))
pbmc <- pbmc[,rownames(meta)]
pbmc@meta.data <- meta
pbmc@reductions <- readRDS(file = paste0("2_pipeline/human_shafflick/", clustering_date, "/reductions_CD4_PBMCs.rds"))
# pbmc <- pbmc[,pbmc$seurat_clusters != 1] ## take out CD8 version integrate_date=2020-04-29
pbmc <- pbmc[,!(pbmc$seurat_clusters %in% c(2, 5, 7, 8))] ## take out non-CD4 version integrate_date=2020-06-01 & clustering_date=2020-06-11

csf <- readRDS(file = paste0("2_pipeline/human_shafflick/", integrate_date, "/so_combined_CSF.rds"))
meta <- readRDS(file = paste0("2_pipeline/human_shafflick/", clustering_date, "/meta_CD4_CSF.rds"))
csf <- csf[,rownames(meta)]
csf@meta.data <- meta
csf@reductions <- readRDS(file = paste0("2_pipeline/human_shafflick/", clustering_date, "/reductions_CD4_CSF.rds"))
csf <- csf[,!(csf$seurat_clusters %in% c(2, 4))]## take out non-CD4 version integrate_date=2020-06-01 & clustering_date=2020-06-11


# #### plot contamination markers
# # load genes of interest and contamination markers
# cont_markers <- read.table("0_data/gene_lists/Genes_for_cleanup.csv", stringsAsFactors = FALSE)[,1]
# cont_markers <- c(cont_markers, "Ptprc", "Tcrg-C4", "Tcrg-C2", "Tcrg-C1") %>% mouse_to_human()
# pdf(file = paste0(dir_out, "UMAP_contamination_pbmc.pdf"), width = 12)
# DefaultAssay(pbmc) <- "RNA"
# pt.size <- 0.5
# for (g in cont_markers) {
#   print(g)
#   if (! (g %in% rownames(pbmc))) {
#     cat(g, ": not in data...\n"); next
#   }
#   print(VlnPlot(object = pbmc, features = g, pt.size = pt.size))
#   print(FeaturePlot(object = pbmc, features = g,
#                     min.cutoff = "q5", max.cutoff = "q95", reduction = "umap", pt.size = pt.size))
# }
# dev.off()
# pdf(file = paste0(dir_out, "UMAP_contamination_csf.pdf"), width = 12)
# DefaultAssay(csf) <- "RNA"
# pt.size <- 0.5
# for (g in cont_markers) {
#   print(g)
#   if (! (g %in% rownames(csf))) {
#     cat(g, ": not in data...\n"); next
#   }
#   print(VlnPlot(object = csf, features = g, pt.size = pt.size))
#   print(FeaturePlot(object = csf, features = g,
#                     min.cutoff = "q5", max.cutoff = "q95", reduction = "umap", pt.size = pt.size))
# }
# dev.off()


DefaultAssay(pbmc) <- "integrated"
DefaultAssay(csf) <- "integrated"
### label possible cluster of th17 cells
pbmc$is_th17 <- pbmc$seurat_clusters %in% c(1)
csf$is_th17 <- csf$seurat_clusters %in% c(0)



##### Signatures ######
## load signature data
bulk_sigs <- readRDS("2_pipeline/spl_clusters/bulk/2020-02-04/signature_slamf6_vs_cxcr6_all.rds")
bulk_sigs <- mouse_to_human(bulk_sigs)
sc_sigs <- readRDS("2_pipeline/differential_expression/SPL_EAE_0_1/2020-04-01/sc_signature.rds")
sc_sigs <- mouse_to_human(sc_sigs)
# from Gaublomme (Meromit)
sigs_raw_gaublomme <- read.table("0_data/gene_lists/all_Gaublomme_2016_signatures.pgf.txt", stringsAsFactors = F)
sigs_raw_gaublomme$sig_name <- paste0(sigs_raw_gaublomme$V1, "-", sigs_raw_gaublomme$V2)
sigs_gaublomme <- lapply(unique(sigs_raw_gaublomme$sig_name), function(i){sigs_raw_gaublomme$V3[sigs_raw_gaublomme$sig_name == i]})
names(sigs_gaublomme) <- unique(sigs_raw_gaublomme$sig_name)
# proliferating from our clustering
sigs_raw_proliferating <- readRDS("2_pipeline/differential_expression/Cluster_UT_GFPpos_inter/2020-04-16/results.rds")
sigs_proliferating <- sigs_raw_proliferating[["7"]] %>%
  rownames_to_column(var = "gene") %>% filter(logFC > log2(1.5) & FDR < 0.05) %>% select(gene) %>% unlist %>% as.character() %>%
  mouse_to_human()

sigs_other <- list()
sigs_other[['Cxcr6_bulk-plus']] <- bulk_sigs$"Cxcr6-plus"
sigs_other[['Slamf6_bulk-plus']] <- bulk_sigs$"Slamf6-plus"
sigs_other[['Cxcr6_sc-plus']] <- sc_sigs$"Cxcr6_sc-plus"
sigs_other[['Slamf6_sc-plus']] <- sc_sigs$"Slamf6_sc-plus"
sigs_other[['Pathogenic_Th17']] <- c("Cxcl3", "Il22", "Il3", "Ccl4", "Gzmb", "Lrmp", "Ccl5", "Casp1", "Csf2", "Ccl3", "Tbx21", "lcos", "ll7r", "Stat4", "Lgals3", "Lag3") %>% mouse_to_human()
sigs_other[["Proliferating"]] <- sigs_proliferating
signature_list <- c(sigs_other, sigs_gaublomme)
saveRDS(signature_list, file = paste0(dir_out, "signature_list.rds"))

# signature_list <- sigs_other

#### compute signature expression score
sapply(signature_list, function(x){
  length(intersect(x, rownames(pbmc)))
})
pbmc <- AddModuleScore(pbmc, features = signature_list, search = F)
colnames(pbmc@meta.data)[grep("Cluster", colnames(pbmc@meta.data))] <- names(signature_list)
pbmc$tissue_MS <- paste(pbmc$tissue, pbmc$is_MS, sep = ".")

sapply(signature_list, function(x){
  length(intersect(x, rownames(csf)))
})
csf <- AddModuleScore(csf, features = signature_list, search = F)
colnames(csf@meta.data)[grep("Cluster", colnames(csf@meta.data))] <- names(signature_list)
csf$tissue_MS <- paste(csf$tissue, csf$is_MS, sep = ".")


saveRDS(pbmc@meta.data, file = paste0(dir_out, "meta_pbmc.rds"))
saveRDS(csf@meta.data, file = paste0(dir_out, "meta_csf.rds"))

#### UMAP all signatures on CD4
pdf(file = paste0(dir_out, "UMAP_signatures_CD4.pdf"))
for (g in names(signature_list)) {
  print(g)
  p1 <- FeaturePlot(pbmc, features = g, min.cutoff = 'q5', max.cutoff = 'q95', split.by = 'tissue_MS', ncol = '2')
  p2 <- FeaturePlot(csf, features = g, min.cutoff = 'q5', max.cutoff = 'q95', split.by = 'tissue_MS', ncol = '2')
  grid.arrange(grobs = list(p1, p2), newpage = TRUE)
}
dev.off()

#### Correlation between Cxcr6 and Th17 pathogenicity
pdf(paste0(dir_out, "correlations_signatures_CD4.pdf"))
plot_df_pbmc <- pbmc@meta.data[,c("tissue_MS", names(signature_list))]
plot_df_csf <- csf@meta.data[,c("tissue_MS", names(signature_list))]

for (g in setdiff(names(signature_list), "Cxcr6_bulk-plus")) {
  print(g)
  plist <- list()
  for (tissue in c("PBMCs", "CSF")) {
    if (tissue == "PBMCs") {
      df_tmp <- plot_df_pbmc[,c("Cxcr6_bulk-plus", g, "tissue_MS")]
    } else if (tissue == "CSF") {
      df_tmp <- plot_df_csf[,c("Cxcr6_bulk-plus", g, "tissue_MS")]
    }
    for (tissue_MS in unique(df_tmp$tissue_MS)) {
      df_tmp_sub <- df_tmp[df_tmp$tissue_MS == tissue_MS,]

      ps_corr <- cor.test(df_tmp_sub[,"Cxcr6_bulk-plus"], df_tmp_sub[,g], alternative = "two.sided", method = "pearson")
      sp_corr <- cor.test(df_tmp_sub[,"Cxcr6_bulk-plus"], df_tmp_sub[,g], alternative = "two.sided", method = "spearman")
      plist[[tissue_MS]] <- ggplot(df_tmp_sub, aes(x = `Cxcr6_bulk-plus`, y = get(g))) +
        geom_point(alpha = 0.6) +
        geom_smooth(color = "blue", method = 'gam') +
        labs(x = "Cxcr6_bulk-plus", y = g) +
        # coord_fixed(ratio = 1) +
        theme_bw() +
        ggtitle(label = tissue_MS, subtitle = paste0("Pearson-corr=", sprintf("%.3f", ps_corr$estimate), ", P-value=", sprintf("%.3f", ps_corr$p.value),
                                                     "\nSpearman-corr=", sprintf("%.3f", sp_corr$estimate), ", P-value=", sprintf("%.3f", sp_corr$p.value)))
    }
  }
  grid.arrange(grobs = plist, nrow = 2)
}
dev.off()



#### Compare Cxcr6 signature score in MS vs. healthy patients in CD4 and Th17

cell_type <- "CD4"
sig_cxcr6 <- c("Cxcr6_bulk-plus", "Cxcr6_sc-plus")

for (cell_type in c("CD4", "Th17")) {
  pdf(file = paste0(dir_out, "Violin_CXCR6_", cell_type, ".pdf"), width = 4, height = 4.5)
  pbmc <- AddModuleScore(pbmc, features = signature_list[sig_cxcr6], search = F, name = "RNA_SIG", assay = "RNA", nbin = 20)
  plot_df <- pbmc@meta.data[, c("orig.ident", "individual", "tissue_MS", "RNA_SIG1", "RNA_SIG2")] %>% rownames_to_column(var = "cell") %>% mutate(tissue_MS = factor(tissue_MS, levels = c("PBMCs.Healthy", "PBMCs.MS")))
  colnames(plot_df)[grep("RNA_SIG", colnames(plot_df))] <- paste0(sig_cxcr6)
  
  for (g in c("Cxcr6_bulk-plus", "Cxcr6_sc-plus")) {
    if (cell_type == "Th17") {
      plot_df <- plot_df[pbmc$is_th17,]
    }
    set.seed(1)
    show_point <- runif(nrow(plot_df)) < 0.2
    p <- ggplot(plot_df, aes(x = individual, y = get(g))) +
      geom_violin(aes(fill = tissue_MS)) +
      geom_jitter(aes(shape = show_point), position=position_jitter(0.4), size = 0.01) +
      scale_shape_manual(values = c(NA, 19), guide = FALSE) +
      scale_fill_manual(values = rev(hcl(h = seq(15, 375, length = nlevels(plot_df$tissue_MS) + 1), l = 65, c = 100)[1:nlevels(plot_df$tissue_MS)])) +
      stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", color = "red", width = 0.5) +
      labs(x = "", y = "Expression Level", fill = "", title = g) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5, face = "bold"),
            panel.background = element_blank(),
            axis.line = element_line()) +
      geom_signif(comparisons = list(c("PBMCs.Healthy", "PBMCs.MS")),
                  step_increase = 0, map_signif_level = T, margin_top = 0.05)
    print(p)
    
    p <- plot_df %>% group_by(individual, tissue_MS) %>% summarise(mean_expr = mean(get(g))) %>% 
      ggplot(aes(x = tissue_MS, y = mean_expr)) +
      # geom_boxplot(aes(fill = tissue_MS)) +
      geom_point(aes(color = tissue_MS)) +
      geom_signif(comparisons = list(c("PBMCs.Healthy", "PBMCs.MS")), test = "t.test") +
      theme_bw()
    print(p)
    
    p <- plot_df %>% 
      ggplot(aes(x = tissue_MS, y = get(g))) +
      geom_boxplot(aes(fill = tissue_MS, group = individual), outlier.shape = NA) +
      ylim(c(-0.5, 0.5)) +
      theme_bw()
    print(p)
    
    p <- plot_df %>% 
      ggplot() +
      geom_density(aes(group = individual, x = get(g), color = tissue_MS), alpha = 0.5) +
      xlim(-0.5, 0.5) +
      theme_bw()
    print(p)
    
    
  }
  dev.off()

}

    
