## configuration
rm(list = setdiff(ls(), "so"))
setwd("/singerlab/linglin/Th17_single_cell_eae_ut/")

library(Seurat)
library(ggplot2)
library(gridExtra)
library(ggsignif)
library(dplyr)
library(tibble)
library(openxlsx)
source("1_code/utils.R")


today <- "2020-07-13"
prep_date <- "2020-06-25"
dir_out <- paste0("2_pipeline/Pdl1_Tibbitt/", today, "/")
if (!dir.exists(dir_out)) {
  dir.create(dir_out, recursive = T)
}


## load data
# so <- readRDS(paste0("2_pipeline/Pdl1_Tibbitt/", prep_date, "/preprocessed_so_RPKM.rds"))

markers <- read.xlsx("0_data/other/Tibbitt_GSE131935_RAW/marker_genes.xlsx", sheet = "6")
isg <- list(markers$gene)
so <- AddModuleScore(so, features = isg, name = "ISG")
VlnPlot(so, features = "ISG1")

DefaultAssay(so) <- "RNA"
so_markers <- FindMarkers(so, ident.1 = 6)
so_markers["Cd274",]
# VlnPlot(so, features = "Cd274", assay = "integrated")
# VlnPlot(so, features = "Cd274", assay = "SCT")
# VlnPlot(so, features = "Cd274", assay = "RNA", split.by = "X.Sample_source_name_ch1", split.plot = F)

# library(edgeR)
# # Create a DGEList data object.
# dgeFull <- DGEList(so[['RNA']]@counts)
# # Estimate the normalization factors
# dgeFull <- calcNormFactors(dgeFull, method="TMM")
# # compute detection rate [here I use scaled number of detected genes]
so@meta.data$cdr <- scale(so$nFeature_RNA)
# # create design matrix
so$in_cluster <- factor(so$seurat_clusters == 6, levels = c(FALSE, TRUE))
# design <- model.matrix(~ in_cluster + X.Sample_source_name_ch1 + cdr, data = droplevels(so@meta.data))
# # Estimate dispersion
# dgeFull <- estimateDisp(dgeFull, design = design)
# # Perform QLF tests
# fit <- glmQLFit(dgeFull, design = design)
# 
# saveRDS(fit$coefficient, file = paste0(dir_out, "edgeR_coefs.rds"))
# saveRDS(fit$dispersion, file = paste0(dir_out, "edgeR_dispersion.rds"))
# 
# 
# qlf <- glmQLFTest(fit, coef = 2)
# res <- topTags(qlf, n = Inf, adjust.method = "BH", sort.by = "logFC", p.value = 1)
# res$table["Cd274",]
# write.csv(res$table, file = paste0(dir_out, "edgeR.csv"))


library(MASS)
tmp <- cbind(so@meta.data[,c("in_cluster", "X.Sample_source_name_ch1")], expr = so[["RNA"]]@data["Cd274",])
tmp %>% group_by(in_cluster, expr > 0) %>% tally() %>% fisher.test()

# glm_fit <- glm(expr ~ in_cluster + X.Sample_source_name_ch1 + cdr, data = tmp,
#                family = negative.binomial(theta = fit$dispersion[which(rownames(fit$coefficient) == "Cd274")], link = "log"))

glm_fit <- glm.nb(expr ~ in_cluster + X.Sample_source_name_ch1 + cdr, data = tmp)
summary(glm_fit)

plot(glm_fit)

t.test(tmp$expr[tmp$in_cluster == TRUE], tmp$expr[tmp$in_cluster == FALSE], paired = F)
tout <- wilcox.test(tmp$expr[tmp$in_cluster == TRUE], tmp$expr[tmp$in_cluster == FALSE], paired = F)
tout$statistic
str(tout)

lm_fit <- lm(expr ~ in_cluster + X.Sample_source_name_ch1 + cdr, data = tmp)
plot(lm_fit)

ggplot(tmp) +
  geom_histogram(aes(x = log(expr + 1), fill = in_cluster), position = "dodge") +
  facet_wrap(~X.Sample_source_name_ch1, nrow = 2)

saveRDS(confint(glm_fit), file = paste0(dir_out, "glm_cd274.rds"))


tmp <- tmp[tmp$X.Sample_source_name_ch1 == "Day 15 BAL T helper cells mouse 2",]
vec <- seq(0, max(tmp$expr), length.out = 100)


tpr <- sapply(vec, function(v) {
  TP <- sum(tmp$in_cluster == TRUE & tmp$expr >= v)
  P <- sum(tmp$in_cluster == TRUE)
  TP/P
})

fpr <- sapply(vec, function(v) {
  FP <- sum(tmp$in_cluster == FALSE & tmp$expr >= v)
  N <- sum(tmp$in_cluster == FALSE)
  FP/N
})

ggplot() +
  geom_path(aes(x = fpr, y = tpr)) +
  # coord_fixed() +
  xlim(0, 1) +
  ylim(0, 1)


# ### create input for COMET
# ## only keep surface protein genes
# surface_genes <- read.table("0_data/gene_lists/COMET_default_genes.txt", sep = ",", header = T) %>% colnames() %>% toupper()
# 
# markers_idx <- match(surface_genes, toupper(rownames(so)))
# markers_idx <- markers_idx[-which(is.na(markers_idx))]
# cat(length(markers_idx), "surface protein genes found!")
# 
# markers_all <- so[["RNA"]]@data[markers_idx,]
# 
# for (tissue in tissue_vec) {
#   dir_out_tissue <- paste0(dir_out, tissue, "/")
#   if (!dir.exists(dir_out_tissue)) dir.create(dir_out_tissue)
#   ### create vis.txt and cluster.txt
#   stopifnot(all(rownames(umap[[tissue]]) == rownames(meta[[tissue]])))
#   ## vis.txt
#   write.table(umap[[tissue]], file = paste0(dir_out_tissue, "vis.txt"), sep = "\t",
#               quote = F, col.names = F, row.names = TRUE)
#   ## cluster.txt
#   cluster_df <- data.frame(
#     cell = rownames(meta[[tissue]]),
#     cluster = meta[[tissue]]$seurat_clusters
#   )
#   write.table(cluster_df, file = paste0(dir_out_tissue, "cluster.txt"), sep = "\t",
#               quote = F, col.names = F, row.names = FALSE)
#   
#   ### create markers.txt
#   markers <- markers_all[,rownames(meta[[tissue]])]
#   write.table(markers, file = paste0(dir_out_tissue, "markers.txt"), sep = "\t",
#               quote = F, row.names = T, col.names = T)
#   
# }