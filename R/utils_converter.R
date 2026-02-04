
symbol2entrez <- function(gene_symbols) {
  all_symbols <- unique(unlist(strsplit(gene_symbols, ";")))
  
  conv_table <- tryCatch({
    bitr(all_symbols, 
         fromType = "SYMBOL", 
         toType = "ENTREZID", 
         OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    return(data.frame(SYMBOL = character(), ENTREZID = character()))
  })
  
  map_genes <- function(x) {
    symbols <- unlist(strsplit(x, ";"))
    entrez <- conv_table$ENTREZID[match(symbols, conv_table$SYMBOL)]
    entrez <- entrez[!is.na(entrez)]
    if(length(entrez) > 0) paste(entrez, collapse = "/") else NA
  }
  
  map_chr(gene_symbols, map_genes)
}

enrichr2enrichResult <- function(df, ontology = "BP", term2gene = NULL) {
  if (is.null(term2gene))
    df <- mutate(df, geneID = symbol2entrez(Genes)) %>% suppressMessages()
  else
    df$geneID <- df$Genes %>% 
      str_replace_all(";", "/") %>% 
      str_split("/") %>% 
      map_chr(~ str_to_sentence(.) %>% paste(collapse = "/"))
  df$Count <- lengths(strsplit(df$geneID, "/"))
  df$ID <- gsub(".*\\((GO:[^)]+)\\).*", "\\1", df$Term)
  df$Term <- gsub(" \\(GO:[^)]+\\)", "", df$Term)
  df$qvalue <- df$Adjusted.P.value
  
  geneSets <- setNames(strsplit(df$geneID, "/"), df$ID)
  genes <- unique(unlist(geneSets))
  gene_universe_size <- length(genes)
  
  df$GeneRatio <- paste0(df$Count, "/", gene_universe_size)
  pathway_ngene <- sapply(strsplit(df$Overlap, "/"), function(x) as.numeric(x[2]))
  if (is.null(term2gene)) {
    db_universe_size <- pathway_ngene %>% sum() %>% divide_by(5) %>% round()
  } else {
    db_universe_size <- pull(term2gene, 2) %>% unique() %>% length()
  }
  df$BgRatio <- paste0(pathway_ngene, "/", db_universe_size)
  
  res <- df[, c("ID", "Term", "GeneRatio", "BgRatio", "P.value", "Adjusted.P.value", "qvalue", "geneID", "Count")]
  colnames(res)[2] <- "Description"
  colnames(res)[5] <- "pvalue"
  colnames(res)[6] <- "p.adjust"
  
  if (is.null(term2gene))
    term2gene <- data.frame(term = rep("dummy", length(genes)), gene = genes)
  else
    term2gene$Pathway <- gsub(" \\(GO:[^)]+\\)", "", term2gene$Pathway)
  
  # if (is.null(clusterProfilerObj))
  clusterProfilerObj <- suppressMessages(
    enricher(
      gene = genes,
      TERM2GENE = term2gene,
      minGSSize = 1
    )
  )
  
  clusterProfilerObj@result <- res
  clusterProfilerObj@gene <- genes
  clusterProfilerObj@geneSets <- geneSets
  clusterProfilerObj@ontology <- ontology
  clusterProfilerObj@organism <- "Homo sapiens"
  clusterProfilerObj@keytype <- "ENTREZID"
  # clusterProfilerObj@pvalueCutoff <- 0.05
  # clusterProfilerObj@pAdjustMethod <- "BH"
  # clusterProfilerObj@method <- character(0) -> clusterProfilerObj@gene2Symbol
  clusterProfilerObj@readable <- FALSE
  # clusterProfilerObj@dr <- list()
  
  return(clusterProfilerObj)
}

enrichr_pairwise_termsim <- function(dummy) {
  geneSets <- dummy@geneSets
  n <- length(geneSets)
  geneNames <- names(geneSets)
  
  
  all_genes <- unique(unlist(geneSets))
  gene_presence <- lapply(geneSets, function(x) all_genes %in% x)
  
  jaccard_mat <- matrix(0, nrow=n, ncol=n)
  rownames(jaccard_mat) <- geneNames
  colnames(jaccard_mat) <- geneNames
  
  for (i in 1:(n-1)) {
    gi <- gene_presence[[i]]
    for (j in (i+1):n) {
      gj <- gene_presence[[j]]
      intersection <- sum(gi & gj)
      union_size <- sum(gi | gj)
      sim <- intersection/union_size
      jaccard_mat[i,j] <- sim
      jaccard_mat[j,i] <- sim
    }
  }
  
  diag(jaccard_mat) <- 1
  
  rownames(jaccard_mat) <- dummy@result$Description -> colnames(jaccard_mat)
  dummy@termsim <- jaccard_mat
  dummy@method <- "JC"
  
  return(dummy)
}

# term2gene_enrichr("CellMarker_2024")
#' @export
term2gene_enrichr <- function(x, to_lower = FALSE) {
  res <- read_delim(
    paste0("https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=", x),
    delim = "\t", 
    col_names = FALSE, 
    trim_ws = TRUE
  )
  
  res <- res %>%
    pivot_longer(
      cols = -X1, 
      names_to = "col", 
      values_to = "Gene"
    ) %>%
    select(Pathway = X1, Gene) %>%
    filter(
      !is.na(Gene),
      Gene != ""
    ) %>%
    mutate(Gene = str_replace(Gene, ",.*$", "")) %>%
    separate_rows(
      Gene,
      sep = "\\s+"
    )
  
  if(to_lower) {
    res <- mutate(res, Gene = str_to_sentence(Gene))
  }
  
  return(res)
}
