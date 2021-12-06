## load and parse GO and KEGG pathways
read_parse_db <- function(f) {
  db_raw <- read.table(f, sep = "@", stringsAsFactors = F)
  db_split <- strsplit(db_raw[,1], split = "\t")
  db <- lapply(db_split, function(x) {
    x[-c(1,2)]
  })
  names(db) <- sapply(db_split, function(x) {x[1]})
  return(db)
}


test_genesets <- function(db, res, reverse = FALSE) {
  ## input for GSEA
  fgesa_stats <- sqrt(res[,'F'])
  fgesa_stats[res$logFC<0] <-  -1 * fgesa_stats[res$logFC<0]
  names(fgesa_stats) <- toupper(rownames(res))
  
  if (reverse) {fgesa_stats = -fgesa_stats}
  
  fgesa_stats <- sort(fgesa_stats)
  
  db <- lapply(db, toupper)
  
  fgsea_0 <- fgsea(pathways = db, 
                   stats = fgesa_stats,
                   nperm = 100000)
  return(fgsea_0)
}
