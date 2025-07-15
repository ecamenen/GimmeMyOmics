bind_gsea <- function(x) {
    x[[1]]@result <- list.map(x, f(i) ~i@result) %>% bind_rows()
    x[[1]]@geneSets <- list.map(x, f(i) ~i@geneSets) %>% flatten()
    return(x[[1]])
}

upset_enrich <- function(
        x,
        regex = NULL,
        width = 20,
        FDR = 0.05,
        ...
) {
    filter_gsea(x, regex = regex, FDR = FDR, width = width) %>%
        upsetplot(...)
}

ridge_gsea <- function(
        x,
        regex = NULL,
        width = 50,
        FDR = 0.05,
        colour = c(palette_discrete()[1], "gray50", palette_discrete()[2]),
        ...
) {
    filter_gsea(x, regex = regex, FDR = FDR, width = width) %>%
        ridgeplot(...)  +
        xlab("NES") %>%
        theme_bulk() %>%
        theme_enrich0(cex = 0.7, colour = colour)
}

#' GSEA enrichment plot with filtering options
#' 
#' Wrapper around `enrichplot::gseaplot2` that adds preprocessing capabilities
#' for GSEA results.
#'
#' @inheritParams enrichplot::gseaplot2
#' @inheritParams network_enrich
#' @param x A `gseaResult` object from clusterProfiler.
#' @param ... Additional arguments passed to `gseaplot2`.
#'
#' @return
#' A ggplot object showing the GSEA enrichment plot for the specified gene set.
#'
#' @examples
#' \dontrun{
#' library(enrichplot)
#' data(geneList, package = "DOSE")
#' gsea_res <- gseGO(geneList, OrgDb = "org.Hs.eg.db", ont = "BP")
#' 
#' # Basic usage (first significant pathway)
#' gsea_enrich(gsea_res)
#' 
#' # With filtering and custom appearance
#' gsea_enrich(
#'   gsea_res,
#'   regex = "signal",
#'   FDR = 0.01,
#'   color = "red",
#'   geneSetID = 2
#' )
#' }
#'
#' @seealso
#' \code{\link[enrichplot]{gseaplot2}}, \code{\link{filter_gsea}}
#'
#' @export
gsea_enrich <- function(
        x,
        regex = NULL,
        width = 20,
        geneSetID = 1,
        color = "black",
        title = NULL,
        FDR = 0.05,
        ...
) {

    filter_gsea(x, regex = regex, FDR = FDR, width = width) %>%
        gseaplot2(
            title = ifelse(is.null(title), .$Description[geneSetID], title),
            geneSetID = geneSetID,
            color = color,
            ...
        )
}

#' Network visualization of enrichment results
#'
#' Wrapper around \code{enrichplot::cnetplot} for enhanced visualization of 
#' gene-concept networks from enrichment analysis results. Provides customization
#' options for node appearance, highlighting specific genes, and improved 
#' aesthetics.
#'
#' @inheritParams enrichplot::cnetplot
#' @param x Either \code{enrichResult} or  \code{compareClusterResult} object from clusterProfiler.
#' @param highlighted Character vector of gene names to highlight in the network.
#' @param col_highlight Color for highlighted nodes.
#' @param cex Double for the scaling factor for all elements.
#' @param cex_node Double for the scaling factor for gene node sizes.
#' @param regex Regular expression pattern to filter by pathway names.
#' @param width Integer for maximum width for pathway name wrapping.
#' @param FDR Double for False Discovery Rate threshold for filtering results.
#' @param colour Color palette for gradient coloring.
#' @param title Character for the plot title.
#' @param power Integer for exponent for color gradient transformation.
#' @param ... Additional arguments passed to \code{enrichplot::cnetplot}.
#'
#' @details
#' This function extends \code{enrichplot::cnetplot} with:
#' \itemize{
#'   \item Automatic filtering by FDR and term description patterns
#'   \item Enhanced label formatting (wrapping, truncation)
#'   \item Customizable node scaling and highlighting
#'   \item Improved color schemes and legends
#' }
#'
#' @examples
#' \dontrun{
#' data(geneList, package = "DOSE")
#' de <- names(geneList)[seq(100)]
#' ora_res <- enrichGO(de, OrgDb = "org.Hs.eg.db", ont = "BP")
#' 
#' # Basic usage
#' network_enrich(ora_res)
#' 
#' # With highlighted genes and custom parameters
#' network_enrich(
#'   ora_res,
#'   foldChange = geneList[seq(100)],
#'   highlighted = c("55143", "991"),
#'   cex = 0.9,
#'   width = 30,
#'   colour = RColorBrewer::brewer.pal(9, "Blues"),
#'   title = "GO Biological Process"
#' )
#' }
#'
#' @return A ggplot object showing the gene-concept network.
#' @seealso \code{\link[enrichplot]{cnetplot}}
#' @export
network_enrich <- function(
        x,
        highlighted = NULL,
        col_highlight = palette_discrete()[1],
        foldChange  = NULL,
        cex = 0.7,
        cex_node = 5.5 * cex,
        node_label = "all",
        regex = NULL,
        width = 20,
        FDR = 0.05,
        colour = brewer.pal(11, "Spectral") %>% rev(),
        title = NULL,
        power = 2,
        ...
) {
    x <- filter_gsea(x, regex = regex, FDR = FDR, width = width)
    # foldChange = rep(1, length(highlighted))
    # names(foldChange) <- highlighted

    if (is(x, "compareClusterResult")) {
        x@compareClusterResult$Description <- to_title(x@compareClusterResult$Description) %>%
            str_trunc(50) %>% str_wrap(width)
    }
    
    p <- cnetplot(
        x,
        node_label = "none",
        cex.params = list(
            category_label = 1.75 * cex,
            gene_label = 1.1 * cex,
            gene_node = cex_node,
            # edge = TRUE,
            category_node = 1.5 * cex
        ),
        color.params = list(foldChange = foldChange, category = "black"),
        shadowtext = "none",
        # colorEdge = TRUE,
        ...
    )
    
    pathway_descriptions <- x@result$Description
    is_pathway <- function(name) name %in% pathway_descriptions
if (!is.null(highlighted)) {
    if (is.null(foldChange)) {
      p <- p + geom_node_text(
        aes_(
          label = ~ifelse(
            .data$name %in% highlighted | is_pathway(.data$name),
            .data$name,
            ""
          ),
          color = ~ifelse(is_pathway(.data$name), "red", "black"),
          size = ~ifelse(is_pathway(.data$name), 1.75 * cex, 0.7 * cex)
        ),
        repel = TRUE,
        show.legend = FALSE
      )
    } else {
      p <- p + geom_node_text(
        aes_(
          label = ~ifelse(
            .data$name %in% highlighted | is_pathway(.data$name),
            .data$name,
            ""
          ),
          size = ~ifelse(is_pathway(.data$name), 1.75 * cex, 0.7 * cex)
        ),
        repel = TRUE,
        show.legend = FALSE
      )
    }
} else {
  
  if (is.null(foldChange)) {
    p <- p + geom_node_text(
      aes_(
        label = ~.data$name,
        color = ~ifelse(is_pathway(.data$name), "red", "black"),
        size = ~ifelse(is_pathway(.data$name), 1.75 * cex, 0.7 * cex)
      ),
      repel = TRUE,
      show.legend = FALSE
    )
  } else {
    p <- p + geom_node_text(
      aes_(
        label = ~.data$name,
        size = ~ifelse(is_pathway(.data$name), 1.75 * cex, 0.7 * cex)
      ),
      repel = TRUE,
      show.legend = FALSE
    )
  }
}
    
    if (is(x, "compareClusterResult")) {
        p <- p +
            scale_fill_manual(
                values = palette_discrete()[seq_along(x@geneClusters)],
                name = ""
            ) +
            guides(color = "none")
    } else {
      title_size <- ifelse(class(x) == "gseaResult", "# Leading genes", "# DEG")
        p <- theme_enrich0(
            p,
            cex * 1.5,
            colour,
            title,
            "Fold change",
            trans = TRUE,
            power = power,
            title_size = title_size,
            num = !is.null(foldChange)
        )
    }

    return(p)
}

heatmap_enrich2 <- function(
    x,
    foldChange,
    wrap = 50,
    regex = NULL,
    width = 20,
    power = 2,
    FDR = 0.05,
    colour = brewer.pal(9, "Spectral") %>% rev(),
    ...
) {
    filter_gsea(x, regex = regex, FDR = FDR, width = width) %>%
        heatplot(
            foldChange = foldChange,
            label_format = function(x) str_trunc(x, wrap) %>% to_title(),
            ...
        ) %>%
        theme_bulk() %>%
        theme_enrich0(
            label_colour = "Fold change",
            colour = colour,
            trans = TRUE,
            power = power
        )
}

tree_enrich <- function(
    x,
    regex = NULL,
    FDR = 0.05,
    metric = "NES",
    width = 20,
    label_words_n = 4,
    group_color = palette_discrete(),
    ...
) {
    filter_gsea(x, regex = regex, FDR = FDR, width = width) %>%
        pairwise_termsim() %>%
        treeplot(
            color = metric,
            cluster.params = list(
                color = group_color[nrow(.) %>% sqrt() %>% seq()],
                label_format = function(x) to_title(x),
                label_words_n = label_words_n
                ),
            ...
            ) %>%
        theme_bulk() +
        scale_colour_gradientn(
            name = metric,
            colours = c(palette_discrete()[1], "gray50", palette_discrete()[2])
        ) +
        guides(
            size = guide_legend(
                order = 1
            ),
            color = guide_colorbar(
                barheight = 5 * cex,
                order = 2,
                frame.colour = "black",
                frame.linewidth = 0.75 * cex,
                ticks.colour = "black",
                ticks.linewidth = 0.75 * cex
            )
        )+
        scale_size_continuous(
            range = c(3, 9),
            breaks = function(x) unique(round(pretty(x, n = 4))),
            labels = label_number(accuracy = 1),
            name = "# Leading genes"
        ) +
        theme(panel.border = element_blank())
}

tree_plot <- function(x, ratio = ratio, wrap = 20, n = 50) {
    # if (is(x, "compareClusterResult")) {
    #   x@compareClusterResult$Description <- to_title(x@compareClusterResult$Description) %>%
    #     str_trunc(wrap)
    # }
    pairwise_termsim(x) %>%
        treeplot(showCategory = 50, offset.params = list(extend = ratio))
}

#' Filter and format GSEA/enrichment results
#'
#' Processes enrichment results by applying FDR cutoff, optional pattern matching,
#' and formatting of term descriptions. Prepares results for visualization.
#'
#' @inheritParams network_enrich
#' @inheritParams stringr::str_detect
#'
#' @details
#' \itemize{
#'   \item Filters results by FDR cutoff (p.adjust ≤ FDR)
#'   \item Optionally filters pathway descriptions using regular expressions
#'   \item Cleans and formats descriptions
#' }
#'
#' @return
#' Returns a filtered and formatted \code{enrichResult} or \code{gseaResult} object 
#' with modified description fields and sorted results.
#'
#' @examples
#' \dontrun{
#' data(geneList, package = "DOSE")
#' 
#' # GSEA example
#' gsea_res <- gseGO(geneList, OrgDb = "org.Hs.eg.db", ont = "BP")
#' filtered_res <- filter_gsea(gsea_res, regex = "signal", FDR = 0.1, width = 30)
#' 
#' # ORA example
#' de_genes <- names(geneList)[1:100]
#' ora_res <- enrichGO(de_genes, OrgDb = "org.Hs.eg.db")
#' filtered_res <- filter_gsea(ora_res, regex = c("apoptosis", "death"), width = 40)
#' }
#'
#' @seealso
#' \code{\link[clusterProfiler]{gseGO}}, \code{\link[clusterProfiler]{enrichGO}}
#'
#' @export
filter_gsea <- function(x, regex = NULL, FDR = 0.05, width = 500, negate = FALSE) {
    x@result <- x@result %>%
        filter(p.adjust <= 0.05) %>%
        mutate(p.adjust = as.numeric(p.adjust))
    if (!is.null(regex)) {
        if (length(regex) > 1) {
            regex <- paste(regex, collapse = "|")
        }
        x@result <- x@result %>%
            filter(str_detect(Description, regex, negate = negate))
    }
    
    new_names <- to_title(x@result[, "Description"]) %>%
      str_remove(" - .*") %>% 
      str_wrap(width)
    
    if (nrow(x@termsim) > 0)
      x@termsim <- x@termsim[x@result[, "Description"], x@result[, "Description"]]
    
    x@geneSets <- x@geneSets[names(x@geneSets) %in% x@result$ID]
    
    x@result[, "Description"] <- new_names
    
    if (nrow(x@termsim) > 0)
      rownames(x@termsim) <- x@result[, "Description"] -> colnames(x@termsim)
    
    if ("NES" %in% colnames(x@result))
      x@result <- x@result %>%
        arrange(NES)
    else 
      x@result <- x@result %>%
        arrange(p.adjust)
    return(x)
}

#' Enrichment map visualization
#' 
#' Wrapper around `enrichplot::emapplot` for creating enrichment maps from 
#' enrichment analysis results. Provides customization options for clustering, 
#' node appearance, and improved aesthetics.
#'
#' @inheritParams enrichplot::emapplot
#' @inheritParams network_enrich
#' @param alpha Double for transparency of edges.
#' @param method A character string specifying the enrichment analysis method. Options are:
#'   - `enrichR` package;
#'   - `clusterProfiler` package.
#' @param ... Additional arguments passed to `enrichplot::emapplot`.
#'
#' @details
#' This function extends `enrichplot::emapplot` with:
#' \itemize{
#'   \item Automatic filtering by FDR and term description patterns
#'   \item Enhanced label formatting (wrapping, truncation)
#'   \item Customizable clustering and scaling parameters
#' }
#'
#' @examples
#' \dontrun{
#' data(geneList, package = "DOSE")
#' de <- names(geneList)[seq(100)]
#' ora_res <- enrichGO(de, OrgDb = "org.Hs.eg.db", ont = "BP")
#' 
#' # Basic usage
#' map_enrich(ora_res)
#' 
#' # With custom parameters
#' map_enrich(
#'   ora_res,
#'   regex = "response",
#'   FDR = 0.01,
#'   cex = 1.2,
#'   showCategory = 20,
#'   colour = c("blue", "white", "red"),
#'   cluster.params = list(cluster = FALSE)
#' )
#' }
#'
#' @return A ggplot object showing the enrichment map.
#' @seealso \code{\link[enrichplot]{emapplot}}, \code{\link[enrichplot]{pairwise_termsim}}
#' @export
map_enrich <- function(
        x,
        regex = NULL,
        width = 20,
        FDR = 0.05,
        cex = 1,
        cluster.params = list(cluster = TRUE, legend = TRUE, label_words_n = 2),
        cex.params = list(category_label =  0.7 * cex, line = 0.5),
        alpha = 0.1,
        showCategory = 15,
        colour = c(palette_discrete()[1], "gray50", palette_discrete()[2]),
        method = "clusterProfiler",
        ...
) {
    p <- filter_gsea(x, regex = regex, FDR = FDR, width = width)
    
    if (nrow(x@termsim) == 0) {
      if (method == "clusterProfiler") {
        p <- pairwise_termsim(p, showCategory = showCategory + 1)
      } else {
        p <- enrichr_pairwise_termsim(p)
      }
    }
    
      p <- emapplot(
          p,
            cluster.params = cluster.params,
            cex.params = cex.params,
            alpha = alpha,
            showCategory = showCategory,
            ...
        ) +
        theme_void() +
        theme(legend.title = element_text(face = "italic", size = 12 * cex))
    
    title_size <- ifelse(class(x) == "gseaResult", "# Leading genes", "# DEG")
    theme_enrich0(p, cex, colour, num = FALSE, title_size = title_size)
}
