library(edgeR)

# edgeR/QLF with cellular detection rate as covariate
run_edgeR_vs_spleen <- function(so) {
  so@meta.data$batch <- factor(so@meta.data$batch, levels = paste0("b", 1:5)) %>% droplevels()
  so@meta.data$tissue <- factor(so@meta.data$tissue, levels = c("SPL", "MLN", "PP", "SI", "COL", "CNS")) %>% droplevels()
  so@meta.data$tissue_batch <- factor(paste(so@meta.data$tissue, so@meta.data$batch, sep = "_"))
  # Create dummy variables
  contrasts(so@meta.data$batch) = contr.sum(nlevels(so@meta.data$batch))
  
  # Create a DGEList data object.
  dgeFull <- DGEList(so@assays$RNA@counts, group=so$tissue, samples = so$orig.ident)
  # Estimate the normalization factors
  dgeFull <- calcNormFactors(dgeFull, method="TMM")
  # compute detection rate [here I use scaled number of detected genes from the original unfiltered data]
  so@meta.data$cdr <- scale(so$nFeature_RNA)
  # create design matrix
  design <- model.matrix(~ tissue + batch + cdr, data = so@meta.data)
  
  # Estimate dispersion
  dgeFull <- estimateDisp(dgeFull, design = design)
  # Perform QLF tests
  fit <- glmQLFit(dgeFull, design = design)
  
  contrast_mat <- matrix(0, nrow = (nlevels(so$tissue_batch) + nlevels(so$tissue) - 1), ncol = ncol(design),
                         dimnames = list(c(levels(so$tissue_batch), setdiff(levels(so$tissue), "SPL")), 
                                         colnames(design)))
  for (tb in rownames(contrast_mat)) {
    tb_vec <- unlist(strsplit(tb, split = "_"))
    if (length(tb_vec) == 2) { # tissue and batch
      if (tb_vec[1] != "SPL") {
        if (tb_vec[2] != "b5") {
          coef_t <- paste0("tissue", tb_vec[1])
          coef_b <- paste0("batch", sub("b", "", tb_vec[2]))
          contrast_mat[tb, c(coef_t, coef_b)] <- 1
        } else {
          coef_t <- paste0("tissue", tb_vec[1])
          contrast_mat[tb, coef_t] <- 1
          coef_b <- paste0("batch", 1:4)
          contrast_mat[tb, coef_b] <- -1
        }
      } else { # other tissues
        if (tb_vec[2] != "b5") {
          coef_b <- paste0("batch", sub("b", "", tb_vec[2]))
          contrast_mat[tb, coef_b] <- 1
        } else {
          coef_b <- paste0("batch", 1:4)
          contrast_mat[tb, coef_b] <- -1
        }
      }
    } else if (length(tb_vec) == 1) { # tissue
      coef_t <- paste0("tissue", tb_vec[1])
      contrast_mat[tb, coef_t] <- 1
    }
  }
  
  
  results <- lapply(rownames(contrast_mat), function(tb){
    cat(tb, "...\n")
    qlf <- glmQLFTest(fit, contrast = contrast_mat[tb,])
    # Extract summaries of differential expression statistics.
    res <- topTags(qlf, n = Inf, adjust.method = "BH", sort.by = "none", p.value = 1)
    return(res$table)
  })
  names(results) <- rownames(contrast_mat)
  results
}




####### process DE results ###########
process_de_vs_spleen_results <- function(so, results, prop_mat, dir_out, thresh_fdr = 0.05, thresh_fc = 1.5){
  # (1) obtain DE genes for each tissue; 
  # (2) compute sample fold changes;
  thresh_lfc <- log2(thresh_fc)
  tissue_vec <- c("MLN", "PP", "SI", "COL")
  # identify DE genes according to specified threshold
  deg <- lapply(tissue_vec, function(y){
    x <- results[[y]] 
    x <- x[order(abs(x$logFC), decreasing = TRUE),] # sort by abs(FC)
    gene <- rownames(x)[x$FDR < thresh_fdr & abs(x$logFC) > thresh_lfc]
    return(data.frame(gene, 
                      "logFC"=x[gene, "logFC"], 
                      "pct_SPL" = prop_mat[gene, "SPL"],  
                      "pct_tissue" = prop_mat[gene, y],
                      "F" = x[gene, "F"],
                      "p_val" = x[gene, "PValue"],
                      "adj_p_val" = x[gene, "FDR"],
                      "logCPM" = x[gene, "logCPM"]))
  })
  names(deg) <- tissue_vec
  deg.up <- lapply(deg, function(x){x[x$logFC>0,]})
  deg.down <- lapply(deg, function(x){x[x$logFC<0,]})
  saveRDS(deg.up, file = paste0(dir_out, "up_all.rds"))
  saveRDS(deg.down, file = paste0(dir_out, "down_all.rds"))
  
  # save DE gene lists
  for (tissue in names(deg)) {
    # up regulated
    deg_table <- deg.up[[tissue]]
    write.table(deg_table,
                file = paste0(dir_out, "up_", tissue, ".csv"),
                row.names = F, sep = ",")
    deg_table_no_rps_rpl  = deg_table %>% filter(!(grepl("^Rps", gene) | grepl("^Rpl", gene)))
    write.table(deg_table_no_rps_rpl,
                file = paste0(dir_out, "up_", tissue, "_no_rps_rpl.csv"),
                row.names = F, sep = ",")
    deg_table_over_70pct  = deg_table %>% filter(pct_tissue > 0.7) # genes that are detected in more than 70% of the cells in the tissue of interest
    write.table(deg_table_over_70pct,
                file = paste0(dir_out, "up_", tissue, "_over_70pct.csv"),
                row.names = F, sep = ",")
    # down regulated
    deg_table <- deg.down[[tissue]]
    write.table(deg_table,
                file = paste0(dir_out, "down_", tissue, ".csv"),
                row.names = F, sep = ",")
    deg_table_no_rps_rpl  = deg_table %>% filter(!(grepl("^Rps", gene) | grepl("^Rpl", gene)))
    write.table(deg_table_no_rps_rpl,
                file = paste0(dir_out, "down_", tissue, "_no_rps_rpl.csv"),
                row.names = F, sep = ",")
    deg_table_les_30pct  = deg_table %>% filter(pct_tissue > 0.3) # genes that are detected in less than 30% of the cells in the tissue of interest
    write.table(deg_table_over_70pct,
                file = paste0(dir_out, "down_", tissue, "_less_30pct.csv"),
                row.names = F, sep = ",")
  }
  
  # compute sample fold changes
  tb_all <- setdiff(names(results), c(tissue_vec, paste0("SPL_b", 1:5)))
  lfc <- sapply(tb_all, function(tb){
    results[[tb]]$logFC
  })
  rownames(lfc) <- rownames(results[[1]])
  saveRDS(lfc, file = paste0(dir_out, "sample_lfc.rds"))
}




# edgeR/QLF with cellular detection rate as covariate
run_edgeR_spl_0_1 <- function(so) {
  so@meta.data$batch <- factor(so@meta.data$batch)
  so@meta.data$seurat_clusters <- factor(so@meta.data$seurat_clusters)
  # Create dummy variables
  contrasts(so@meta.data$batch) = contr.sum(nlevels(so@meta.data$batch))
  
  # Create a DGEList data object.
  dgeFull <- DGEList(so[['RNA']]@counts)
  # Estimate the normalization factors
  dgeFull <- calcNormFactors(dgeFull, method="TMM")
  # compute detection rate [here I use scaled number of detected genes]
  so@meta.data$cdr <- scale(so$nFeature_RNA)
  # create design matrix
  design <- model.matrix(~ seurat_clusters + batch + cdr, data = so@meta.data)
  
  # Estimate dispersion
  dgeFull <- estimateDisp(dgeFull, design = design)
  # Perform QLF tests
  fit <- glmQLFit(dgeFull, design = design)
  qlf <- glmQLFTest(fit, coef = 2)
  res <- topTags(qlf, n = Inf, adjust.method = "BH", sort.by = "none", p.value = 1)
  return(res)
}


### each cluster vs. all others
run_edgeR_clusters <- function(so) {
  so@meta.data$batch <- factor(so@meta.data$batch)
  so@meta.data$seurat_clusters <- factor(so@meta.data$seurat_clusters)
  # Create dummy variables
  contrasts(so@meta.data$batch) = contr.sum(nlevels(so@meta.data$batch))
  
  # Create a DGEList data object.
  dgeFull <- DGEList(so[['RNA']]@counts)
  # Estimate the normalization factors
  dgeFull <- calcNormFactors(dgeFull, method="TMM")
  # compute detection rate [here I use scaled number of detected genes]
  so@meta.data$cdr <- scale(so$nFeature_RNA)
  
  res <- list()
  for (cl in levels(so$seurat_clusters)) {
    cat("DE for cluster", cl, "\n")
    so$in_cluster <- so$seurat_clusters == cl
    # create design matrix
    design <- model.matrix(~ in_cluster + batch + cdr, data = so@meta.data)
    # Estimate dispersion
    dgeFull <- estimateDisp(dgeFull, design = design)
    # Perform QLF tests
    fit <- glmQLFit(dgeFull, design = design)
    qlf <- glmQLFTest(fit, coef = 2)
    res[[cl]] <- topTags(qlf, n = Inf, adjust.method = "BH", sort.by = "none", p.value = 1)$table
  }
  
  return(res)
}





run_edgeR_vs_spleen_tissue_trt <- function(so) {
  so$tissue_treatment <- factor(so$tissue_treatment) %>% relevel(ref = "SPL_UT")
  results <- list()
  for (tissue_treatment in levels(so$tissue_treatment)) {
    if (tissue_treatment == "SPL_UT") next()
    if (tissue_treatment == "DLN_EAE") next()
    cat(tissue_treatment, "\n")
    so_sub <- so[,so$tissue_treatment %in% c("SPL_UT", tissue_treatment)]
    so_sub$batch <- factor(so_sub$batch) %>% droplevels()
    so_sub$tissue_treatment <- droplevels(so_sub$tissue_treatment)
    # Create dummy variables
    contrasts(so_sub$batch) = contr.sum(nlevels(so_sub$batch))
    # Create a DGEList data object.
    dgeFull <- DGEList(so_sub@assays$RNA@counts)
    # Estimate the normalization factors
    dgeFull <- calcNormFactors(dgeFull, method="TMM")
    # compute detection rate [here I use scaled number of detected genes from the original unfiltered data]
    so_sub$cdr <- scale(so_sub$nFeature_RNA)
    # create design matrix
    design <- model.matrix(~ tissue_treatment + batch + cdr, data = so_sub@meta.data)
    
    # Estimate dispersion
    dgeFull <- estimateDisp(dgeFull, design = design)
    # Perform QLF tests
    fit <- glmQLFit(dgeFull, design = design)
    # saveRDS(fit, file = paste0(dir_out, "glmQLFit.rds"))
    qlf <- glmQLFTest(fit, coef = 2)
    # Extract summaries of differential expression statistics.
    res <- topTags(qlf, n = Inf, adjust.method = "BH", sort.by = "none", p.value = 1)
    saveRDS(res, file = paste0(dir_out, tissue_treatment, "_results.rds"))
    results[[tissue_treatment]] <- res$table
    
  }
  return(results)
}



# # edgeR/QLF with cellular detection rate as covariate
# run_edgeR_tfh_gfppos_vs_gfpneg <- function(so) {
#   so$batch <- factor(so$batch)
#   so$GFP_positive <- factor(so$GFP_positive, levels = c(FALSE, TRUE))
#   # Create dummy variables
#   contrasts(so$batch) = contr.sum(nlevels(so$batch))
#   
#   # Create a DGEList data object.
#   dgeFull <- DGEList(so[['RNA']]@counts)
#   # Estimate the normalization factors
#   dgeFull <- calcNormFactors(dgeFull, method="TMM")
#   # compute detection rate [here I use scaled number of detected genes]
#   so@meta.data$cdr <- scale(so$nFeature_RNA)
#   # create design matrix
#   design <- model.matrix(~ GFP_positive + batch + cdr, data = so@meta.data)
#   # Estimate dispersion
#   dgeFull <- estimateDisp(dgeFull, design = design)
#   # Perform QLF tests
#   fit <- glmQLFit(dgeFull, design = design)
#   qlf <- glmQLFTest(fit, coef = 2)
#   res <- topTags(qlf, n = Inf, adjust.method = "BH", sort.by = "none", p.value = 1)
#   return(res)
# }


