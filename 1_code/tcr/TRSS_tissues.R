rm(list = setdiff(ls(), "so_all"))

setwd("/singerlab/linglin/Th17_single_cell_eae_ut")
library(ggplot2)
library(RColorBrewer)
library(Seurat)
library(dplyr)
library(tibble)
library(tidyr)
library(entropy)
library(grid)
library(gridExtra)
library(ggpubr)
library(corrplot)
library(pheatmap)
library(openxlsx)
library(dendextend)

set.seed(1)
source("1_code/utils.R")
source("1_code/tcr/utils.R")

#### configuration ####
cargs <- commandArgs(trailingOnly = TRUE)
if (length(cargs) == 0) {
  today <- Sys.Date()#"2020-12-10"
  clonotyping_type_date <- "dominant_2020-04-02"
  prep_type_date <- "dominant_TCR_2020-03-24"
} else {
  today <- cargs[1]
  clonotyping_type_date <- cargs[2]
  prep_type_date <- cargs[3]
}

dir_out <- paste0("2_pipeline/TCR/compare_tissue/")
if(!dir.exists(dir_out)) dir.create(dir_out, recursive = T)

#### load data ####
# tissue_vec <- c("SPL", "MLN", "PP", "SI", "COL", "CNS", "DLN")

## meta data
meta_all <- readRDS(paste0("for_paper/results/preprocessing//meta_processed_", prep_type_date, ".rds"))

## clonotype data
clty <- read.csv(paste0("2_pipeline/TCR/clonotype_assignment_", clonotyping_type_date, ".csv"), row.names = 1)

## merge data
meta_all$clonotype_id <- clty[rownames(meta_all), "clonotype_id"]
meta_all$pooled_clonotype_id <- paste(meta_all$treatment, meta_all$batch, meta_all$mouse, meta_all$clonotype_id, sep = "_")
meta_noNA <- meta_all %>% filter(!is.na(clonotype_id))

#### compare tissues ###
treatment_vec <- c("EAE", "UT")
batch_vec <- c("b6", "b7", "b8", "b9")
mouse_vec <- c("m1", "m2")

for (treatment in treatment_vec) {
  cl_mtx <- list()
  ## pooled
  cl_mtx[["pooled"]] <- with(meta_noNA[meta_noNA$treatment == treatment,], table(pooled_clonotype_id, tissue) %>% as.matrix())
  ## each mouse
  for (batch in batch_vec) {
    for (mouse in mouse_vec) {
      meta_noNA_sub <- meta_noNA[meta_noNA$treatment == treatment & meta_noNA$batch == batch & meta_noNA$mouse == mouse, ]
      if (nrow(meta_noNA_sub) == 0) next ## mouse does not exist
      sample_name <- paste(treatment, batch, mouse, sep = "_")
      cl_mtx[[sample_name]] <- with(meta_noNA_sub, table(clonotype_id, tissue) %>% as.matrix())
    }
  }

  cl_mtx <- lapply(cl_mtx, function(x){
    apply(x, 2, function(y){y/sum(y)})
  })
  bh_dist <- lapply(cl_mtx, function(x){
    tmp <- philentropy::distance(t(x), method = "bhattacharyya") %>% as.matrix()
    rownames(tmp) <- colnames(x)
    colnames(tmp) <- colnames(x)
    tmp
  })

  bh_dist_cap <- bh_dist
  bh_dist_cap <- lapply(bh_dist_cap, function(x){
    pmin(x, 10*max(x[!is.infinite(x)]))
  })
  hc <- lapply(bh_dist_cap, function(x){hclust(as.dist(x), method = "ward.D2")})
  saveRDS(hc, paste0(dir_out, "bh_hclust_", treatment, "_", today, ".rds"))

  bh_sim <- lapply(bh_dist, function(x){exp(-x)})
  bh_sim_plt <- lapply(names(bh_sim), function(i){
    x <- bh_sim[[i]]
    diag(x) <- NA
    x[hc[[i]]$order, hc[[i]]$order]
  })
  names(bh_sim_plt) <- names(bh_sim)
  
  saveRDS(bh_sim_plt, paste0(dir_out, "bh_similarity_", treatment, "_", today, ".rds"))

  # pdf(paste0(dir_out, "similarity_", treatment, ".pdf"))
  # for (i in names(bh_sim_plt)) {
  #   layout(matrix(c(1, 1, 2, 2, 2, 2, 2, 2, 2, 2), nrow = 10, ncol = 1, byrow = TRUE))
  #   
  #   ## dendrogram
  #   par(mar=c(0, 7, 2, 6.5))
  #   # par(fig = c(0,10,6.5,10)/10)
  #   dend <- as.dendrogram(hc[[i]])
  #   dend %>% set("labels", "") %>%
  #     plot(main = paste0(treatment, " ", i))
  # 
  #   ## similarity
  #   # par(fig = c(0,10,0,6.5)/10)
  #   # par(new = TRUE)
  #   corrplot(bh_sim_plt[[i]],
  #            is.corr = F,
  #            type="upper", order="original", method = "circle",
  #            # col=rev(brewer.pal(n=8, name="RdYlBu")),
  #            outline = T,
  #            col = rev(col_corr(200)),
  #            mar = c(0, 0, 0, 0),
  #            tl.col = tissue_colors(rownames(bh_sim_plt[[i]])),
  #            tl.pos = 'td', tl.cex = 2,
  #            diag = T, na.label = " ")
  # }
  # dev.off()
}



#### heatmap by clonotypes ####
meta_noNA$sample_name <- paste(meta_noNA$treatment, meta_noNA$batch, meta_noNA$mouse, sep = "_") %>% factor()
sample_name_vec <- c("pooled", levels(meta_noNA$sample_name))

for (treatment in treatment_vec) {
  plt_mtx_list <- list()
  for (i in sample_name_vec) {
    if (i == "pooled") {
      tmp <- meta_noNA[meta_noNA$treatment == treatment,]
      tmp$clonotype_id <- tmp$pooled_clonotype_id
    } else {
      tmp <- meta_noNA[meta_noNA$treatment == treatment & meta_noNA$sample_name == i,]
      if (nrow(tmp) == 0) next()
    }
    cl_to_plot <- tmp %>% group_by(clonotype_id) %>% tally() %>% arrange(-n)
    plt_mtx <- table(tmp$clonotype_id, tmp$tissue) %>% as.matrix()
    plt_mtx <- plt_mtx[cl_to_plot$clonotype_id,]
    plt_mtx <- plt_mtx[rowSums(plt_mtx > 0) > 0,]
    plt_mtx_list[[i]] <- plt_mtx
  }
  saveRDS(plt_mtx_list, paste0(dir_out, "size_table_", treatment, "_", today, ".rds"))
}


# min_size_vec <- c(1, 5, 10)
# meta_noNA$sample_name <- paste(meta_noNA$treatment, meta_noNA$batch, meta_noNA$mouse, sep = "_") %>% factor()
# sample_name_vec <- c("pooled", levels(meta_noNA$sample_name))
# plt_scale_type_vec <- c("col_nocap", "col_cap", "none_cap") #c("col_cap", "row_nocap", "none_cap")
# for (min_size in min_size_vec) {
#   for (treatment in treatment_vec) {
#     pdf(paste0(dir_out, "heatmap_clonotype_", treatment, "_min", min_size, ".pdf"), width = 4, height = 6)
#     for (i in sample_name_vec) {
#       if (i == "pooled") {
#         tmp <- meta_noNA[meta_noNA$treatment == treatment,]
#         tmp$clonotype_id <- tmp$pooled_clonotype_id
#       } else {
#         tmp <- meta_noNA[meta_noNA$treatment == treatment & meta_noNA$sample_name == i,]
#         if (nrow(tmp) == 0) next()
#       }
#       cl_to_plot <- tmp %>% group_by(clonotype_id) %>% tally() %>% arrange(-n) %>% filter(n >= min_size)
#       plt_mtx <- table(tmp$clonotype_id, tmp$tissue) %>% as.matrix()
#       plt_mtx <- plt_mtx[cl_to_plot$clonotype_id,]
#       plt_mtx <- plt_mtx[rowSums(plt_mtx > 0) > 0,]
#       for (plt_scale_type in plt_scale_type_vec) {
#         if (plt_scale_type == "col_nocap") {
#           plt_mtx_scale <- 100 * apply(plt_mtx, 2, function(x){x/sum(x)})
#           plt_title <- paste0(i, " (column scaled as %, not capped)")
#           cluster_rows <- TRUE
#           cluster_cols <- TRUE
#         } else if (plt_scale_type == "col_cap") {
#           plt_mtx_scale <- 100 * apply(plt_mtx, 2, function(x){x/sum(x)})
#           ceiling_value <- quantile(plt_mtx_scale, 0.95)
#           plt_mtx_scale <- pmin(plt_mtx_scale, ceiling_value)
#           plt_title <- paste0(i, " (column scaled as %, capped at ", round(ceiling_value, 2), ")")
#           cluster_rows <- TRUE
#           cluster_cols <- TRUE
#         } else if (plt_scale_type == "row_nocap") {
#           plt_mtx_scale <- 100 * apply(plt_mtx, 1, function(x){x/sum(x)}) %>% t
#           plt_title <- paste0(i, " (row scaled as %, not capped)")
#           cluster_rows <- TRUE
#           cluster_cols <- TRUE
#         } else if (plt_scale_type == "none_cap") {
#           plt_mtx_scale <- plt_mtx
#           ceiling_value <- quantile(plt_mtx_scale, 0.99)
#           plt_mtx_scale <- pmin(plt_mtx_scale, ceiling_value)
#           # n_mouse <- rowSums(plt_mtx>0)
#           # n_cell <- rowSums(plt_mtx)
#           # plt_mtx_scale <- plt_mtx_scale[order(n_mouse, n_cell, decreasing = T),]
#           plt_title <- paste0(i, " (not scaled, capped at ", round(ceiling_value, 2), ")")
#           cluster_rows <- TRUE
#           cluster_cols <- TRUE
#         }
#         
#         p <- pheatmap(plt_mtx_scale, 
#                       cluster_rows = cluster_rows, 
#                       cluster_cols = cluster_cols, 
#                       clustering_method = "ward.D2",
#                       show_rownames = F, 
#                       color = rev(col_corr_pos_skip(100)),
#                       main = plt_title, 
#                       fontsize = 8,
#                       fontsize_col = 11,
#                       silent = T)
#         which_col_names <- which(p$gtable$layout$name == "col_names")
#         p$gtable$grobs[[which_col_names]]$gp$col = tissue_colors(p$tree_col$labels[p$tree_col$order])
#         p$gtable$grobs[[which_col_names]]$rot = 0
#         p$gtable$grobs[[which_col_names]]$hjust = 0.5
#         p$gtable$grobs[[which_col_names]]$vjust = 1
#         grid.arrange(p$gtable)
#       }
#     }
#     dev.off()
#     
#   }
# }



  
