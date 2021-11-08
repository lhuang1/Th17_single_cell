######## plot overview after clustering #######
plot_cluster_overview <- function(so, dir_out, tissue = NULL, pt.size = 0.5, 
                                  pie_by = NULL, 
                                  cluster_cols = NULL,
                                  cols_by = c("seurat_clusters", "batch")){
  if (is.null(tissue)) {
    pdf(file = paste0(dir_out, "/PLOTS/UMAP_overview.pdf"))
  } else {
    pdf(file = paste0(dir_out, "/PLOTS/", tissue, "/UMAP_overview.pdf"))
  }
  
  # sample sizes
  p <- so@meta.data %>%
    group_by(orig.ident) %>% 
    tally() %>% 
    ggplot(aes(x=orig.ident, y = n))+
    geom_bar(aes(fill = orig.ident), stat = "identity") +
    geom_text(aes(y = n, label = n)) +
    ggtitle(paste0("Sample sizes (Total = ", ncol(so), " cells)")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
  
  # cluster sizes
  p <- so@meta.data %>%
    group_by(seurat_clusters) %>% 
    tally() %>%
    ggplot(aes(x=seurat_clusters, y = n))+
    geom_bar(aes(fill = seurat_clusters), stat = "identity") +
    geom_text(aes(y = n, label = n), vjust=0) +
    theme_bw() +
    ggtitle(paste0("Cluster sizes (Total = ", ncol(so), " cells)"))
  print(p)
  
  # UMAP by given factors
  for (x in cols_by) {
    if (x %in% colnames(so@meta.data)) {
      x_ = x
    } else {
      cat(x, "not found in meta data...\n")
      next
    }
    
    p <- DimPlot(object = so, reduction = "umap", group.by = x_, pt.size = pt.size, label = FALSE)
    p$data <- p$data[sample(1:nrow(p$data), nrow(p$data), replace = F),]
    
    if (x == "GFP_positive") {
      p <- p + scale_color_manual(values = c("lightgrey", "darkgreen"), 
                                  limits = c(FALSE, TRUE),
                                  name = "", labels = c('Ex Th17', 'Current Th17'))
    } else if (x == "tissue") {
      p <- p + scale_color_manual(values = tissue_colors(unique(so$tissue)))
    } else if (x == "treatment") {
      p$layers[[1]]$aes_params$alpha <- 0.7
    } 
    print(p)
  }
  
  # pie charts for cluster compositions
  if (!is.null(pie_by)){
    par(mfrow = c(3,2))
    for (pct_by in pie_by) {
      plist <- plot_pie_charts(so, pct_by)
      for (p in plist) {print(p)}
    }
  }
  
  dev.off()
}


### Make pie charts
plot_pie_charts <- function(so, pct_by){
  slices_all <- table(so@meta.data[, pct_by])
  if (pct_by == "GFP_positive") {
    col_vec <- colorRampPalette(c("lightgrey", "darkgreen"))(2)
  } else if (pct_by == "tissue") {
    col_vec <- tissue_colors(unique(so$tissue))
  } else {
    col_vec <- rainbow(length(slices_all))
  }
  
  cluster_names <- unique(so@active.ident)
  plist <- list()
  for (x in cluster_names) {
    slices_x <- table(so@meta.data[so@active.ident == x, pct_by])
    slices <- slices_x[names(slices_all)]
    slices[is.na(slices)] <- 0 ## if does not exist in current cluster, set to 0
    pct <- round(slices/sum(slices)*100)
    lbls <- paste0(names(slices), "=", pct, "%") # add percents to labels 
    plist[[paste0(x, "_unscaled")]] <- pie(slices, labels = lbls, col = col_vec,
                      main=paste0("Cluster ", x, ", total cells: ", sum(slices), ", unscaled"))
    
    slices_scaled <- slices /slices_all
    pct_scaled <- round(slices_scaled/sum(slices_scaled)*100)
    lbls_scaled <- paste0(names(slices_scaled), "=", pct_scaled, "%") # add percents to labels 
    plist[[paste0(x, "_scaled")]] <- pie(slices_scaled, labels = lbls_scaled, col = col_vec,
                      main=paste0("Cluster ", x, ", total cells: ", sum(slices), ", scaled"))
  }
  return(plist)
}


save_and_plot_markers <- function(so, so_markers, dir_markers, tissue = NULL, pt.size = 0.1){
  if (is.null(tissue)) {
    dir_markers_file <- paste0(dir_markers, "/FILES/")
    dir_markers_plot <- paste0(dir_markers, "/PLOTS/")
  } else {
    dir_markers_file <- paste0(dir_markers, "/FILES/", tissue, "/")
    dir_markers_plot <- paste0(dir_markers, "/PLOTS/", tissue, "/")
  }
  if (!dir.exists(dir_markers_file)) dir.create(dir_markers_file, recursive = T)
  if (!dir.exists(dir_markers_plot)) dir.create(dir_markers_plot, recursive = T)
  
  cluster_names <- unique(so_markers$cluster)
  for (x in cluster_names) {
    cat("Cluster", x, "...\n")
    cluster_markers = so_markers %>% filter(cluster == x) #%>% arrange(desc(avg_logFC))
    write.table(cluster_markers, 
                file = paste0(dir_markers_file, "Cluster_Markers_", x, ".txt"), quote = F)
    cluster_markers_no_rps_rpl  = cluster_markers %>% 
      filter(!(grepl("^Rps", gene) | grepl("^Rpl", gene)))
    write.table(cluster_markers_no_rps_rpl, 
                file = paste0(dir_markers_file, "Cluster_Markers_", x, "_no_rps_rpl.txt"), quote = F)
    cluster_markers_thresh03  = cluster_markers %>% 
      filter(pct.2 < 0.3)
    write.table(cluster_markers_thresh03, 
                file = paste0(dir_markers_file, "Cluster_Markers_", x, "_thresh03.txt"), quote = F)
    to_plot <- cluster_markers$gene[1:min(18, nrow(cluster_markers))]
    pdf(file = paste0(dir_markers_plot, "Cluster_Markers_", so@active.assay, "_", x, ".pdf"))
    for (to_plot_i in to_plot) {
      # to_plot_i <- ((i-1)*6 + 1) : (i*6)
      print(VlnPlot(object = so, features = to_plot_i, pt.size = pt.size))
      print(FeaturePlot(object = so, features = to_plot_i, 
                        min.cutoff = "q5", 
                        max.cutoff = "q95",
                        reduction = "umap", pt.size = pt.size))
    }
    dev.off()
  }
  return(NULL)
}
