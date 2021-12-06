gene_lists_to_modules <- function(deg_list, collapse_mLN = TRUE){
  genes <- Reduce(union, lapply(deg_list, function(x){x$gene}))
  de_mat <- sapply(deg_list, function(x){genes %in% x$gene})
  rownames(de_mat) <- genes
  n_de_vec <- rowSums(de_mat)
  stopifnot(all(n_de_vec>0))
  de_combn <- lapply(ncol(de_mat):1, function(n_de){combn(names(deg_list), n_de)})
  modules <- lapply(de_combn, function(x){
    tmp <- lapply(1:ncol(x), function(y){
      z <- x[, y]
      if (length(z) == 1){
        rownames(de_mat)[n_de_vec == length(z) & as.integer(de_mat[,z]) == length(z)]
      } else {
        rownames(de_mat)[n_de_vec == length(z) & rowSums(de_mat[,z]) == length(z)]
      }
    })
    names(tmp) <- apply(x, 2, function(xx){paste(xx, collapse = "_")})
    return(tmp)
  })
  
  if (collapse_mLN) { # only keep "mLN only" group; collapse other combinations involving mLN to reduce categories
    modules[[4]]$PP <- union(modules[[4]]$PP, modules[[3]]$mLN_PP); modules[[3]]$mLN_PP <- NULL
    modules[[4]]$SI <- union(modules[[4]]$SI, modules[[3]]$mLN_SI); modules[[3]]$mLN_SI <- NULL
    modules[[4]]$colon <- union(modules[[4]]$colon, modules[[3]]$mLN_colon); modules[[3]]$mLN_colon <- NULL
    modules[[3]]$PP_SI <- union(modules[[3]]$PP_SI, modules[[2]]$mLN_PP_SI); modules[[2]]$mLN_PP_SI <- NULL
    modules[[3]]$PP_colon <- union(modules[[3]]$PP_colon, modules[[2]]$mLN_PP_colon); modules[[2]]$mLN_PP_colon <- NULL
    modules[[3]]$SI_colon <- union(modules[[3]]$SI_colon, modules[[2]]$mLN_SI_colon); modules[[2]]$mLN_SI_colon <- NULL
    modules[[2]]$PP_SI_colon <- union(modules[[2]]$PP_SI_colon, modules[[1]]$mLN_PP_SI_colon); modules[[1]] <- NULL
  }
  
  return(modules)
}
