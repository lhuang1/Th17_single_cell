## mapping between hashtags (assigned my demuxEM) to tissue label for batch 6, 7 and 8
hashtag2tissue_b6_b7_b8 <- function(idents){
  c("301_302" = "SPL",
    "303_304" = "MLN",
    "301_305" = "PP",
    "302_305" = "SI",
    "302_303" = "COL",
    "301_304" = "CNS",
    "304_305" = "DLN"
  )[idents]
}
## mapping between hashtags (assigned my demuxEM) to tissue label for batch 9
hashtag2tissue_b9 <- function(idents){
  c("301_302" = "SPL",
    "303_304" = "MLN",
    "301_305" = "PP",
    "303_305" = "SI",
    "302_304" = "COL"
  )[idents]
}

## parse demuxEM assignment (hashtags label) to tissue labels
parse_assignment <- function(ident, batch) {
  library(dplyr)
  library(stringr)
  tissue_label <- rep(NA, length(ident))
  n_tags <- str_count(ident, pattern = ",")
  tissue_label[n_tags == 0] <- "Negative"
  tissue_label[n_tags >= 2] <- "Multiplet" 
  idx_valid <- which(n_tags == 1)
  tags_valid <- strsplit(ident[idx_valid], split = ",") %>% 
    sapply(., function(x){sort(as.numeric(x)) %>% paste(collapse = "_")})
  if (batch == "b9") {
    tissue_label[idx_valid] <- hashtag2tissue_b9(tags_valid)
  } else {
    tissue_label[idx_valid] <- hashtag2tissue_b6_b7_b8(tags_valid)
  }
  tissue_label[is.na(tissue_label)] <- "Other"
  
  all_labels <- c("SPL", "MLN", "PP", "SI", "COL", "CNS", "DLN", "Negative", "Multiplet", "Other")
  tissue_label <- factor(tissue_label, levels = all_labels) %>% droplevels()
  
  return(tissue_label)
}


### plot cell QC measures for a Seurat object (so) ####
plot_cell_QC <- function(so){
  p <- VlnPlot(object = so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  print(p)
  p <- FeatureScatter(object = so, feature1 = "nCount_RNA", feature2 = "percent.mt")
  print(p)
  p <- FeatureScatter(object = so, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(p)
}

#### main function for expression data preprocessing ####
prep_data <- function(tissue, treatment, batch, mouse, project_name, dir_in){
  # pdf(paste0(dir_out, "QC_", project_name, ".pdf"))
  
  if (tissue != "hashed") { # non-hashing
    count_dir <- paste0(dir_in, batch, "_", tissue, "_", treatment, "_filtered_feature_bc_matrix/")
    gfp_dir <- paste0(dir_in, batch, "_", tissue, "_", treatment, "_GFP_raw_feature_bc_matrix/")
  } else { # hashing
    count_dir <- paste0(dir_in, batch, "_hashed_", treatment, "_", mouse, "_filtered_feature_bc_matrix/")
    gfp_dir <- paste0(dir_in, batch, "_hashed_", treatment, "_", mouse, "_GFP_raw_feature_bc_matrix/")
  }
  
  
  if (!dir.exists(count_dir)) return(NULL)
  if (!dir.exists(gfp_dir)) {
    print("GFP data not available!!!!")
    stop()
  }
  
  ## load gene counts and combine with GFP
  cts <- Read10X(count_dir, strip.suffix = T)
  gfp_cts <- Read10X(gfp_dir, strip.suffix = T)
  rownames(gfp_cts)[rownames(gfp_cts) == "IL17A_GFP_wg"] <- "Il17a-GFP"
  cts <- rbind(cts, "Il17a-GFP" = gfp_cts["Il17a-GFP", colnames(cts)])
  rm(gfp_cts) # free memory
  
  so <- CreateSeuratObject(cts, project = project_name, min.cells = 0, min.features = 1)
  
  #### cell QC
  ncell_before <- ncol(so)
  so[["percent.mt"]] <- PercentageFeatureSet(object = so, pattern = "^mt-") # compute mitochondrial gene percentage
  plot_cell_QC(so) # QC plots before filtering out cells
  
  # compute proportion of house keeping genes that are detected in each cell
  so[["HK_prop"]] <- Matrix::colMeans(so[["RNA"]]@counts[intersect(hk_genes, rownames(so)),] > 0)
  p <- ggplot() +
    geom_histogram(aes(x = so$HK_prop), color = "white") +
    scale_y_log10()# histogram should be bimodal
  print(p)
  
  # add tissue info
  if (tissue == "hashed") {
    meta_dir <- paste0(dir_proj, "output/results/preprocessing/demultiplex/meta_", batch, "_hashed_", treatment, "_", mouse, ".rds")
    meta <- readRDS(meta_dir) 
    if (treatment == "UT") {
      meta <- meta %>% filter(assignment_tissue %in% c("SPL", "PP", "MLN", "SI", "COL")) %>% droplevels()
    } else {
      meta <- meta %>% filter(assignment_tissue %in% c("SPL", "PP", "MLN", "SI", "COL", "CNS", "DLN")) %>% droplevels()
    }
    so$tissue <- meta$assignment_tissue[match(gsub("-.*", "", colnames(so)), rownames(meta))]
    so <- so[,!is.na(so$tissue)]
  } else {
    so$tissue <- tissue
  }
  
  ## subset data
  so <- subset(x = so, subset = nFeature_RNA > 500 & percent.mt < 5 & HK_prop > 0.6)
  
  ncell_after <- ncol(so)
  plot_cell_QC(so) # QC plots after filtering out cells
  
  if (ncell_after < 100) return(NULL) # do not return if <100 cells passed QC
  
  # plot number of cells before and after cell QC
  tmp <- data.frame(
    stage =  factor(c("before", "after"), levels = c("before", "after")),
    cell_number = c(ncell_before, ncell_after)
  )
  ggplot(tmp, aes(x = stage, y = cell_number)) + 
    geom_bar(stat = "identity", width=0.5) +
    geom_text(aes(label = cell_number), vjust = 1.6, color="white", size = 3.5)
  # dev.off()
  
  ## add other meta info
  so$batch <- batch
  so$treatment <- treatment
  so$is_hashing <- tissue == "hashed"
  so$GFP_positive <- so[['RNA']]@counts["Il17a-GFP",] > 0
  
  if (tissue == "hashed") {
    ## split by tissue
    so_by_tissue <- SplitObject(so, split.by = "tissue")
    so_by_tissue <- lapply(so_by_tissue, function(x){
      x$orig.ident <- paste(x$tissue, x$orig.ident, sep = "_")
      x
    })
    names(so_by_tissue) <- paste0(names(so_by_tissue), "_", project_name)
    so_by_tissue[sapply(so_by_tissue, ncol) < 50] <- NULL # do not return if <50 cells in a tissue passed QC
  } else {
    so_by_tissue <- list(so)
    names(so_by_tissue) <- project_name
  }
  
  return(so_by_tissue)
}


####### Identify cluster(s) that are specific to one batch #########
## Algorithm overview:
# For each cluster, 
# a. find the 2 batch with highest in cluster fraction (n in cluster / n in and not in)
# b. Test if top 1 is significantly higher than top 2 with chi square 
# c. If significant, check odds ratio. 
#   i. OR > 2, mark as batch-specific
#   ii. OR <= 2, not batch-specific
# d. If not significant, not batch
find_batch_cluster <- function(meta) {
  bc <- c() ## batch clusters in current sample
  bsize_all <- meta %>% group_by(batch) %>% tally()
  
  cluster_vec <- sort(unique(meta$seurat_clusters))
  for (cl in cluster_vec) {
    meta_sub <- meta[meta$seurat_clusters == cl,]
    bsize_in <- meta_sub %>% group_by(batch) %>% tally()
    
    if (nrow(bsize_in) == 1) {
      bc <- c(bc, cl) ## must be a batch specific cluster if only one batch found
      next
    }
    
    bsize_df <- merge(bsize_all, bsize_in, by = 1, all = TRUE)
    bsize_df[is.na(bsize_df)] <- 0 ## replace NAs with 0s
    bsize_df$n.z <- bsize_df$n.x - bsize_df$n.y
    bfrac_in <- bsize_df$n.y/bsize_df$n.x
    top_2 <- order(bfrac_in, decreasing = T)[1:2]
    
    t_out <- fisher.test(bsize_df[top_2, c("n.y", "n.z")]) # Fisher's exact test
    if (t_out$p.value < 0.05 &  abs(log2(t_out$estimate)) > 2) bc <- c(bc, cl)
  }
  return(bc)
}

    



