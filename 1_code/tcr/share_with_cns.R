rm(list = setdiff(ls(), "so_all"))

setwd("/singerlab/linglin/Th17_single_cell_eae_ut")
library(ggplot2)
library(RColorBrewer)
library(Seurat)
library(dplyr)
library(tibble)
library(tidyr)
library(ggpubr)
library(edgeR)
library(pheatmap)
library(gridExtra)
set.seed(1)
source("1_code/utils.R")
source("1_code/tcr/utils.R")
source("1_code/misc/egghatch.R")

#### configuration ####
cargs <- commandArgs(trailingOnly = TRUE)
if (length(cargs) == 0) {
  today <- "2020-08-03"
  clonotyping_type_date <- "dominant_2020-04-02"
  prep_type_date <- "dominant_TCR_2020-03-24"
} else {
  today <- cargs[1]
  clonotyping_type_date <- cargs[2]
  prep_type_date <- cargs[3]
}

dir_out <- paste0("2_pipeline/TCR/share_with_CNS/", today, "/")
if(!dir.exists(dir_out)) dir.create(dir_out, recursive = T)

#### load data ####
tissue_vec <- c("MLN", "PP", "SI", "COL", "SPL", "CNS", "DLN")
tissue_vec_sub <- c("MLN", "PP", "SI", "COL")
# so_all <- readRDS(paste0("2_pipeline/preprocessing/so_processed_", prep_type_date, ".rds"))
# so_all@reductions <- readRDS(paste0("2_pipeline/clustering/all_inter/2020-04-16/FILES/reductions.rds"))
so <- so_all[,so_all$treatment == "EAE"]

## reductions
reductions_ls <- lapply(tissue_vec, function(tissue) {
  if (tissue == "SPL") {
    reductions <- readRDS(paste0("2_pipeline/clustering/EAE_GFPall_SPL/2020-03-25/FILES/reductions.rds"))
  } else if (tissue == "DLN") {
    reductions <- readRDS(paste0("2_pipeline/clustering/all_intra/2020-05-07/FILES/DLN/reductions.rds"))
  } else {
    reductions <- readRDS(paste0("2_pipeline/clustering/EAE_GFPall_intra/2020-03-25/FILES/", tissue, "/reductions.rds"))
  }
})
names(reductions_ls) <- tissue_vec


## clonotype data
clty <- read.csv(paste0("2_pipeline/TCR/clonotype_assignment_", clonotyping_type_date, ".csv"), row.names = 1)

## merge data
so$clonotype_id <- clty[colnames(so), "clonotype_id"]
so$pooled_clonotype_id <- paste(so$treatment, so$batch, so$mouse, so$clonotype_id, sep = "_")
meta_noNA <- so@meta.data %>% rownames_to_column(var = "cell") %>% filter(!is.na(clonotype_id))

## find clones shared with CNS
meta_noNA <- meta_noNA %>% group_by(pooled_clonotype_id) %>% 
  mutate(shared_with_CNS = any(tissue == "CNS"))

cells_shared_with_CNS <- sapply(tissue_vec, function(tissue) {
  meta_noNA$cell[meta_noNA$tissue == tissue & meta_noNA$shared_with_CNS] %>% as.character()
})


### plot number of cells sharing
plt_df <- meta_noNA %>% group_by(tissue) %>% 
  summarise(n_shared = sum(shared_with_CNS), frac_shared = mean(shared_with_CNS))

p_nshared <- plt_df %>% filter(tissue %in% tissue_vec_sub) %>% 
  mutate(tissue = factor(tissue, levels = tissue_vec_sub)) %>% 
ggplot() +
  geom_bar(aes(x = tissue, y =  n_shared, fill = tissue), stat = "identity", width = 0.8, show.legend = FALSE) +
  geom_text(aes(x = tissue, y = n_shared/2, label = paste0(n_shared, "\n(", sprintf("%.1f", 100*frac_shared), "%)"))) +
  scale_fill_manual(values = tissue_colors()) +
  labs(x = "", y = "Number of Cells") +
  theme_bw()

pdf(paste0(dir_out, "n_cells_share_with_CNS.pdf"), width = 3, height = 3)
p_nshared
dev.off()


### check sharing status with SPL_0 vs. SPL_1
meta_spl <- read.table("2_pipeline/clustering/EAE_GFPall_SPL/2020-03-25/FILES/meta_data.txt", sep = " ")
cl_spl_0 <- meta_noNA$pooled_clonotype_id[match(rownames(meta_spl)[meta_spl$seurat_clusters == 0 & !is.na(meta_spl$raw_clonotype_id)], meta_noNA$cell)] %>% as.character()
cl_spl_1 <- meta_noNA$pooled_clonotype_id[match(rownames(meta_spl)[meta_spl$seurat_clusters == 1 & !is.na(meta_spl$raw_clonotype_id)], meta_noNA$cell)] %>% as.character()
cl_spl_0_1 <- intersect(cl_spl_0, cl_spl_1)
cl_spl_all <- meta_noNA$pooled_clonotype_id[match(rownames(meta_spl)[!is.na(meta_spl$raw_clonotype_id)], meta_noNA$cell)] %>% as.character()

## for cells that share TCR with CNS, which SPL cluster do they share clonotype with?
plt_df <- sapply(cells_shared_with_CNS[tissue_vec_sub], function(x){
  in_spl_0 <- meta_noNA$pooled_clonotype_id[match(x, meta_noNA$cell)] %in% setdiff(cl_spl_0, cl_spl_0_1)
  in_spl_1 <- meta_noNA$pooled_clonotype_id[match(x, meta_noNA$cell)] %in% setdiff(cl_spl_1, cl_spl_0_1)
  in_spl_0_1 <- meta_noNA$pooled_clonotype_id[match(x, meta_noNA$cell)] %in% cl_spl_0_1
  in_spl_other <- meta_noNA$pooled_clonotype_id[match(x, meta_noNA$cell)] %in% setdiff(cl_spl_all, union(cl_spl_0, cl_spl_1))
  not_in_spl <- !(meta_noNA$pooled_clonotype_id[match(x, meta_noNA$cell)] %in% cl_spl_all)
  c(n_in_spl_0 = sum(in_spl_0), 
    # frac_in_spl_0 = mean(in_spl_0), 
    n_in_spl_1 = sum(in_spl_1), 
    # frac_in_spl_1 = mean(in_spl_1),
    n_in_spl_0_1 = sum(in_spl_0_1), 
    n_in_spl_other = sum(in_spl_other),
    n_not_in_spl = sum(not_in_spl)
    )
}) %>% t %>% data.frame() %>% tibble::rownames_to_column(var = "tissue") %>% 
  gather(key = "type", value = "n_cells", -tissue) %>% 
  mutate(tissue = factor(tissue, levels = tissue_vec)) %>% 
  mutate(type = factor(type, 
                       levels = c("n_in_spl_0", "n_in_spl_1", "n_in_spl_0_1", "n_in_spl_other", "n_not_in_spl"), 
                       labels = c("SPL_0", "SPL_1", "SPL_0&1", "SPL_other", "Not found"))) 

p_cell_n <- plt_df %>% 
  # filter(type %in% c("n_in_spl_0", "n_in_spl_1")) %>% 
  ggplot() +
  geom_bar(aes(x = tissue, y = n_cells, fill = type), stat = "identity") +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  labs(y = "Number of cells", x = "", fill = "")


p_cell_pct <- plt_df %>% 
  group_by(tissue) %>% mutate(pct = 100 * n_cells / sum(n_cells)) %>%
  # filter(type %in% c("n_in_spl_0", "n_in_spl_1")) %>% 
  ggplot() +
  geom_bar(aes(x = tissue, y = pct, fill = type), stat = "identity") +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  labs(y = "Percent of cells", x = "", fill = "")



## for clonotypes shared with CNS, in which SPL cluster were they also found?
plt_df <- sapply(cells_shared_with_CNS[tissue_vec_sub], function(x){
  in_spl_0 <- intersect(meta_noNA$pooled_clonotype_id[match(x, meta_noNA$cell)], setdiff(cl_spl_0, cl_spl_0_1))
  in_spl_1 <- intersect(meta_noNA$pooled_clonotype_id[match(x, meta_noNA$cell)], setdiff(cl_spl_1, cl_spl_0_1))
  in_spl_0_1 <- intersect(meta_noNA$pooled_clonotype_id[match(x, meta_noNA$cell)], cl_spl_0_1)
  in_spl_other <- intersect(meta_noNA$pooled_clonotype_id[match(x, meta_noNA$cell)], setdiff(cl_spl_all, union(cl_spl_0, cl_spl_1)))
  not_in_spl <- setdiff(meta_noNA$pooled_clonotype_id[match(x, meta_noNA$cell)], cl_spl_all)
  
  c(n_in_spl_0 = length(in_spl_0), 
    n_in_spl_1 = length(in_spl_1), 
    n_in_spl_0_1 = length(in_spl_0_1), 
    n_in_spl_other = length(in_spl_other),
    n_not_in_spl = length(not_in_spl)
  )
}) %>% t %>% data.frame() %>% tibble::rownames_to_column(var = "tissue")

p_clone_n <- plt_df %>% gather(key = "type", value = "n", -tissue) %>% 
  mutate(tissue = factor(tissue, levels = tissue_vec)) %>% 
  mutate(type = factor(type, 
                       levels = c("n_in_spl_0", "n_in_spl_1", "n_in_spl_0_1", "n_in_spl_other", "n_not_in_spl"), 
                       labels = c("SPL_0", "SPL_1", "SPL_0&1", "SPL_other", "Not found"))) %>% 
  ggplot() +
  geom_bar(aes(x = tissue, y = n, fill = type), stat = "identity") +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  labs(y = "Number of clonotypes", x = "", fill = "")
p_clone_pct <- plt_df %>% gather(key = "type", value = "n", -tissue) %>% 
  mutate(tissue = factor(tissue, levels = tissue_vec)) %>% 
  group_by(tissue) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup() %>% 
  mutate(type = factor(type, 
                       levels = c("n_in_spl_0", "n_in_spl_1", "n_in_spl_0_1", "n_in_spl_other", "n_not_in_spl"), 
                       labels = c("SPL_0", "SPL_1", "SPL_0&1", "SPL_other", "Not found"))) %>% 
  ggplot() +
  geom_bar(aes(x = tissue, y = pct, fill = type), stat = "identity") +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  labs(y = "Percent of clonotypes", x = "", fill = "")



## expansion of clones in SPL that are shared by tissue, SPL and CNS
plt_df <- sapply(cells_shared_with_CNS[tissue_vec_sub], function(x){
  in_spl_0 <- intersect(meta_noNA$pooled_clonotype_id[match(x, meta_noNA$cell)], setdiff(cl_spl_0, cl_spl_0_1))
  in_spl_1 <- intersect(meta_noNA$pooled_clonotype_id[match(x, meta_noNA$cell)], setdiff(cl_spl_1, cl_spl_0_1))
  in_spl_0_1 <- intersect(meta_noNA$pooled_clonotype_id[match(x, meta_noNA$cell)], cl_spl_0_1)
  # in_spl_other <- intersect(meta_noNA$pooled_clonotype_id[match(x, meta_noNA$cell)], setdiff(cl_spl_all, union(cl_spl_0, cl_spl_1)))
  # not_in_spl <- setdiff(meta_noNA$pooled_clonotype_id[match(x, meta_noNA$cell)], cl_spl_all)
  tmp_0 <- meta_noNA %>% 
    filter(tissue == "SPL" & (pooled_clonotype_id %in% in_spl_0)) %>% 
    group_by(pooled_clonotype_id) %>% tally()
  
  tmp_1 <- meta_noNA %>% 
    filter(tissue == "SPL" & (pooled_clonotype_id %in% in_spl_1)) %>% 
    group_by(pooled_clonotype_id) %>% tally()

  tmp_0_1 <- meta_noNA %>% 
    filter(tissue == "SPL" & (pooled_clonotype_id %in% in_spl_0_1)) %>% 
    group_by(pooled_clonotype_id) %>% tally()  
  
  
  c("SPL_0" = mean(tmp_0$n), "SPL_1" = mean(tmp_1$n), "SPL_0_1" = mean(tmp_0_1$n))
  
})

p_exp_shared <- pheatmap(plt_df, 
                         color = colorRampPalette(colors = c("white", "red"))(10),
                         breaks = seq(0, 20, length.out = 11),
                         cluster_rows = F, 
                         cluster_cols = F,
                         display_numbers = T,
                         number_color = "black",
                         fontsize_number = 20,
                         main = "Average clone size share by\neach tissue & CNS & SPL")


pdf(paste0(dir_out, "share_with_SPL_clusters_", today, ".pdf"), width = 4, height = 4)
p_cell_n
p_cell_pct
p_clone_n
p_clone_pct
grid.arrange(p_exp_shared$gtable)
dev.off()


### plot umap
plist <- list()
for (tissue in tissue_vec) {
  if (tissue == "CNS") next()
  umap_df <- reductions_ls[[tissue]]$umap@cell.embeddings %>% data.frame()
  umap_df$shared_with_CNS <- rownames(umap_df) %in% cells_shared_with_CNS[[tissue]]
  umap_df$shared_with_CNS[!(rownames(umap_df) %in% meta_noNA$cell)] <- NA
  plist[[tissue]] <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(aes(color = shared_with_CNS), size = 0.5) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
    ggtitle(tissue) +
    labs(x = "", y = "") +
    theme_classic() +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0.5))
}

tmp <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = shared_with_CNS), size = 1) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), na.value = "grey60") +
  theme_classic()
p_legend <- get_legend(tmp)

pdf(paste0(dir_out, "UMAP_shared.pdf"), width = 7, height = 9)
grid.arrange(grobs = plist, right = p_legend, left = text_grob("UMAP_2", rot = 90, size = 10), bottom = text_grob("UMAP_1", rot = 0, size = 10))
dev.off()



### plot PC
plist <- list()
for (tissue in tissue_vec) {
  if (tissue == "CNS") next()
  pca_df <- reductions_ls[[tissue]]$pca@cell.embeddings %>% data.frame()
  pca_df$shared_with_CNS <- rownames(pca_df) %in% cells_shared_with_CNS[[tissue]]
  pca_df$shared_with_CNS[!(rownames(pca_df) %in% meta_noNA$cell)] <- NA
  plist[[tissue]] <- ggplot(pca_df, aes(x = PC_1, y = PC_2)) +
    geom_point(aes(color = shared_with_CNS), size = 0.5, alpha = 0.6) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), na.value = "grey60") +
    ggtitle(tissue) +
    labs(x = "", y = "") +
    theme_classic() +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0.5))
}

tmp <- ggplot(pca_df, aes(x = PC_1, y = PC_2)) +
  geom_point(aes(color = shared_with_CNS), size = 1) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  theme_classic()
p_legend <- get_legend(tmp)

pdf(paste0(dir_out, "pca_shared.pdf"), width = 7, height = 9)
grid.arrange(grobs = plist, right = p_legend, left = text_grob("PC_2", rot = 90, size = 10), bottom = text_grob("PC_1", rot = 0, size = 10))
dev.off()




# ### differential expression
# min_cells <- 500
# res <- list()
# for (tissue in tissue_vec) {
#   if (tissue == "CNS") next()
#   cat(tissue, "...\n")
#   
#   so_sub <- so[,so$tissue == tissue]
#   so_sub$shared_with_CNS <- meta_noNA$shared_with_CNS[match(colnames(so_sub), meta_noNA$cell)]
#   so_sub <- FindNeighbors(so_sub, assay = "integrated", k.param = 50)
#   # so_sub <- FindClusters(so_sub, resolution = 2)
#   # so_sub <- RunUMAP(so_sub, dims = 1:15)
#   # DimPlot(so_sub)
#   # DimPlot(so_sub, group.by = "shared_with_CNS")
#   
#   not_shared <- meta_noNA$cell[meta_noNA$tissue == tissue & !meta_noNA$shared_with_CNS]
#   if (length(cells_shared_with_CNS[[tissue]]) < min_cells) {# not shared less than min_cells
#     add_not_shared <- which(Matrix::colSums(so_sub@graphs$integrated_nn[not_shared,]) > 0) %>% names()
#     add_shared <- which(Matrix::colSums(so_sub@graphs$integrated_nn[cells_shared_with_CNS[[tissue]],]) > 0) %>% names()
#   } else {
#     add_not_shared <- c()
#     add_shared <- c()
#   }
#   
#   so_sub@meta.data[setdiff(add_shared, add_not_shared), "shared_with_CNS"] <- TRUE
#   # so_sub@meta.data[setdiff(add_not_shared, add_shared), "shared_with_CNS"] <- FALSE
#   meta_sub <- so_sub@meta.data[!is.na(so_sub$shared_with_CNS),]
#   meta_sub$batch_mouse <- factor(paste(meta_sub$batch, meta_sub$mouse, sep = "_"))
#   
#   cts_mtx <- so[["RNA"]]@counts[,rownames(meta_sub)]
#   frac_0 <- Matrix::rowMeans(cts_mtx == 0)
#   cts_mtx <- cts_mtx[frac_0 < 0.9,]
# 
#   # Create a DGEList data object.
#   dgeFull <- DGEList(cts_mtx)
#   # Estimate the normalization factors
#   dgeFull <- calcNormFactors(dgeFull, method="TMM")
#   # compute detection rate [here I use scaled number of detected genes from the original unfiltered data]
#   meta_sub$cdr <- scale(meta_sub$nFeature_RNA)
#   # create design matrix
#   design <- model.matrix(~ shared_with_CNS + batch_mouse + cdr, data = meta_sub)
#   # Estimate dispersion
#   cat("\testimating dispersion...\n")
#   dgeFull <- estimateDisp(dgeFull, design = design)
#   # Perform QLF tests
#   fit <- glmQLFit(dgeFull, design = design)
#   qlf <- glmQLFTest(fit, coef = 2)
#   res[[tissue]] <- topTags(qlf, n = Inf, adjust.method = "BH", sort.by = "none", p.value = 1)
# }
# saveRDS(res, paste0(dir_out, "res_edgeR.rds"))
# 
# # res <- readRDS(paste0(dir_out, "res_edgeR.rds"))
# up_genes <- lapply(res, function(x){
#   rownames(x)[x$table$FDR < 0.05 & x$table$logFC > log2(1.5)]
# })
# all_up_genes <- unique(unlist(up_genes))
# up_mtx <- sapply(names(res), function(tissue){
#   all_up_genes %in% up_genes[[tissue]]
# })
# rownames(up_mtx) <- all_up_genes
# 
# pdf(paste0(dir_out, "venn_deg.pdf"))
# vennDiagram(vennCounts(up_mtx[,-6]))
# dev.off()


### check PC space
### plot PC
pc_max <- list()
plist <- list()
pca_df <- list()
pdf(paste0(dir_out, "pca_shared_top2.pdf"), width = 8, height = 9)
for (tissue in tissue_vec_sub) {
  pca_df[[tissue]] <- reductions_ls[[tissue]]$pca@cell.embeddings %>% data.frame()
  pca_df[[tissue]]$shared_with_CNS <- rownames(pca_df[[tissue]]) %in% cells_shared_with_CNS[[tissue]]
  pca_df[[tissue]]$shared_with_CNS[!(rownames(pca_df[[tissue]]) %in% meta_noNA$cell)] <- NA
  
  library(pROC)
  pca_df_sub <- pca_df[[tissue]] %>% filter(complete.cases(pca_df[[tissue]]))
  auc_vec <- sapply(colnames(reductions_ls[[tissue]]$pca@cell.embeddings), function(i) {
    rocobj <- roc(pca_df_sub$shared_with_CNS, pca_df_sub[,i])
    abs(as.numeric(rocobj$auc) - 0.5)
  }) %>% sort(decreasing = T)
  
  pc_max[[tissue]] = names(auc_vec)[1:2]
  
  plist[[tissue]] <-
  ggplot(pca_df[[tissue]], aes(x = get(pc_max[[tissue]][1]), y = get(pc_max[[tissue]][2]))) +
    geom_point(aes(color = shared_with_CNS), size = 0.5, alpha = 0.6) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), na.value = "grey60") +
    ggtitle(tissue) +
    labs(x = pc_max[[tissue]][1], y = pc_max[[tissue]][2]) +
    theme_classic() +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0.5))
  print(plist[[tissue]])
  
}
dev.off()

## MLN, pc_2 +
## SI, pc_6 +
## PP, pc_1 +
mln <- names(sort(reductions_ls[["MLN"]]$pca@feature.loadings[,"PC_2"]))[1:200]
si <- names(sort(reductions_ls[["SI"]]$pca@feature.loadings[,"PC_6"]))[1:200]
pp <- names(sort(reductions_ls[["PP"]]$pca@feature.loadings[,"PC_1"]))[1:200]

vc_mtx <- sapply(list(mln, si, pp), function(x){
  Reduce(union, list(mln,si,pp)) %in% x
})
colnames(vc_mtx) <- c("MLN", "SI", "PP")

vennDiagram(vennCounts(vc_mtx))


## use signature
res <- readRDS(paste0("2_pipeline/TCR/share_with_CNS/2020-07-19/res_edgeR.rds"))
up <- res$SI$table %>% filter(FDR < 0.05 & logFC > log2(1.5)) %>% 
  rownames_to_column(var = "gene") %>%
  arrange(PValue)
write.csv(up, file = paste0(dir_out, "SI_CNS_sharing_signature_", today, ".csv"))

sigs <- list(SI_up = up$gene)

# volcano plot
library(ggrepel)
p_volcano <- res$SI$table %>% 
  rownames_to_column(var = "gene") %>% 
  mutate(is_de = FDR < 0.05 & abs(logFC) > log2(1.5)) %>% 
  mutate(gene_label = ifelse(PValue < 1e-10 & abs(logFC) > 1, gene, NA)) %>% 
  ggplot(aes(x = logFC, y = -log10(PValue), color = is_de)) +
  geom_point() +
  geom_label_repel(aes(label = gene_label), show.legend = FALSE) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
  theme_bw()
pdf(paste0(dir_out, "SI_volcano_", today, ".pdf"), width = 10, height = 7)
p_volcano
dev.off()


meta_sub_ls <- list()
for (tissue in tissue_vec_sub) {
  so_sub <- so[, rownames(reductions_ls[[tissue]]$pca@cell.embeddings)]
  so_sub <- AddModuleScore(so_sub, features = sigs, name = "SI_CNS")
  so_sub$shared_with_CNS <- meta_noNA$shared_with_CNS[match(colnames(so_sub), meta_noNA$cell)]
  so_sub$shared_with_CNS[is.na(so_sub$shared_with_CNS)] <- "No TCR"
  so_sub$shared_with_CNS <- factor(so_sub$shared_with_CNS, levels = c("TRUE", "FALSE", "No TCR"), labels = c("Shared", "Not shared", "No TCR"))
  meta_sub_ls[[tissue]] <- so_sub@meta.data
}

sig_scores <- Reduce(rbind, meta_sub_ls) %>% 
  rownames_to_column(var = "cell_name") %>% 
  select(cell_name, score = SI_CNS1, tissue, shared_with_CNS)
saveRDS(sig_scores, paste0(dir_out, "SI_CNS_sharing_signature_score_", today, ".rds"))

library(ggsignif)

p_score_vln <- sig_scores %>% 
  filter(shared_with_CNS %in% c("Shared", "Not shared")) %>% 
  ggplot(aes(x = shared_with_CNS, y = score)) +
  geom_violin(aes(fill = tissue)) +
  facet_wrap(~tissue, nrow = 1) +
  geom_signif(test="t.test", comparisons = list(c("Shared", "Not shared")), map_signif_level = T) +
  scale_fill_manual(values = tissue_colors()) +
  theme_bw()


p_score_box <- sig_scores %>% 
  filter(shared_with_CNS %in% c("Shared", "Not shared")) %>% 
  ggplot(aes(x = shared_with_CNS, y = score)) +
  geom_boxplot(aes(fill = tissue), show.legend = F, width = 0.5) +
  facet_wrap(~tissue, nrow = 1) +
  geom_signif(test="t.test", comparisons = list(c("Shared", "Not shared")), map_signif_level = T) +
  scale_fill_manual(values = tissue_colors()) +
  labs(x = "", y = "SI-CNS sharing signature score", fill = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(paste0(dir_out, "SI_CNS_sharing_signature_score_", today, ".pdf"), width = 5, height = 5)
p_score_vln
p_score_box
dev.off()


#### gsea for SI shared vs. not shared
## DB
source("1_code/enrichment/utils.R")
library(fgsea)
go_bp <- read_parse_db("0_data/other/GO_Biological_Process_2018.txt")
go_bp <- human_to_mouse(go_bp) %>% lapply(., toupper)

kegg <- read_parse_db("0_data/other/KEGG_2019_Mouse.txt")

## input for GSEA
gsea <- list()
gsea$go_bp <- test_genesets(db = go_bp, res = res$SI$table, reverse = FALSE)
gsea$kegg <- test_genesets(db = kegg, res = res$SI$table, reverse = FALSE)

library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "go_bp")
writeData(wb, sheet = "go_bp", gsea$go_bp[order(gsea$go_bp$padj),], rowNames = F, colNames = T)
addWorksheet(wb, "kegg")
writeData(wb, sheet = "kegg", gsea$kegg[order(gsea$kegg$padj),], rowNames = F, colNames = T)
saveWorkbook(wb, file = paste0(dir_out, "GSEA_GO_KEGG.xlsx"), overwrite = T)

