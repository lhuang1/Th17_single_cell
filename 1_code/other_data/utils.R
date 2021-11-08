save_and_plot_markers <- function(so, so_markers, dir_markers, tissue = NULL, pt.size = 0.1){
  if (is.null(tissue)) {
    dir_markers_file <- paste0(dir_markers, "/FILES/")
    dir_markers_plot <- paste0(dir_markers, "/PLOTS/")
  } else {
    dir_markers_file <- paste0(dir_markers, "/FILES/", tissue, "/")
    dir_markers_plot <- paste0(dir_markers, "/PLOTS/", tissue, "/")
  }
  dir.create(dir_markers_file, recursive = T)
  dir.create(dir_markers_plot, recursive = T)
  
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