bind_gsea <- function(x) {
    x[[1]]@result <- list.map(x, f(i) ~i@result) %>% bind_rows()
    x[[1]]@geneSets <- list.map(x, f(i) ~i@geneSets) %>% flatten()
    return(x[[1]])
}

upset_enrich <- function(
        x,
        regex = NULL,
        width = 20,
        ...
) {
    filter_gsea(x, regex = regex, FDR = FDR, width = width) %>%
        upsetplot(...)
}

ridge_enrich <- function(
        x,
        regex = NULL,
        width = 20,
        colour = c(GimmeMyPlot:::palette_discrete()[1], "gray50", GimmeMyPlot:::palette_discrete()[2]),
        ...
) {
    filter_gsea(x, regex = regex, FDR = FDR, width = width) %>%
        ridgeplot(...) %>%
        theme_bulk() %>%
        theme_enrich0(cex = 0.7, colour = colour) +
        xlab("NES")
}

gsea_enrich <- function(
        x,
        regex = NULL,
        width = 20,
        geneSetID = 1,
        id = 1,
        color = "black",
        title = NULL,
        ...
) {

    filter_gsea(x, regex = regex, FDR = FDR, width = width) %>%
        gseaplot2(title = ifelse(is.null(title), .$Description[id], title), geneSetID = geneSetID, color = color, ...)
}

network_enrich <- function(
        x,
        highlighted = NULL,
        col_highlight = GimmeMyPlot:::palette_discrete()[1],
        foldChange  = NULL,
        cex = 0.7,
        cex_node = 5.5 * cex,
        node_label = "all",
        regex = NULL,
        width = 20,
        FDR = 0.05,
        # colour = c(rep("grey", 6), "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142"),
        colour = brewer.pal(11, "Spectral") %>% rev(),
        title = NULL,
        power = 2,
        title_size = "# Leading genes",
        ...
) {
    x <- filter_gsea(x, regex = regex, FDR = FDR, width = width)
    # foldChange = rep(1, length(highlighted))
    # names(foldChange) <- highlighted

    if (is(x, "compareClusterResult")) {
        x@compareClusterResult$Description <- to_title(x@compareClusterResult$Description) %>%
            str_trunc(50) %>% str_wrap(width)
    }

    # options(ggrepel.max.overlaps = 1000)
    p <- cnetplot(
        x,
        node_label = node_label,
        cex.params = list(
            category_label = 1.75 * cex,
            gene_label = 1.1 * cex,
            gene_node = cex_node,
            # edge = TRUE,
            category_node = 1.5 * cex
        ),
        # max.overlaps = Inf,
        color.params = list(foldChange = foldChange, category = "black"),
        shadowtext = "none",
        # colorEdge = TRUE,
        ...
    )
    if (is(x, "compareClusterResult")) {
        p <- p +
            scale_fill_manual(
                values = GimmeMyPlot:::palette_discrete()[seq_along(x@geneClusters)],
                name = ""
            ) +
            guides(color = "none")
    } else {
        p <- theme_enrich0(
            p,
            cex * 1.5,
            colour,
            title,
            "Fold change",
            trans = TRUE,
            power = power,
            title_size = title_size,
            num = TRUE
        )
        # scale_color_gradientn(
        #   # trans =
        # ) +
    }

    if (!is.null(highlighted)) {
        # pathways <- is.na(p$data$color)
        highlighted0 <- p$data$name %in% highlighted
        p$data$size[highlighted0] <- 50
        # p$data$name[!highlighted0 & ! pathways] <- ""
        p$data$name[!highlighted0] <- ""
    }

    return(p)
}

heatmap_enrich2 <- function(x, foldChange, wrap = 50, regex = NUL, width = 20, power = 2, colour = brewer.pal(9, "Spectral") %>% rev(), ...) {
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


tree_enrich <- function(x, regex = NULL, FDR = 0.05, metric = "NES", width = 20, label_words_n = 4, group_color = GimmeMyPlot:::palette_discrete(), ...) {
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
            colours = c(GimmeMyPlot:::palette_discrete()[1], "gray50", GimmeMyPlot:::palette_discrete()[2])
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
            labels = scales::label_number(accuracy = 1),
            name = "# Leading genes"
        ) +
        theme(panel.border = element_blank())
}

tree_plot <- function(x, ratio = ratio, wrap = 20, n = 50) {
    # if (is(x, "compareClusterResult")) {
    #   x@compareClusterResult$Description <- to_title(x@compareClusterResult$Description) %>%
    #     str_trunc(wrap)
    # }
    enrichplot::pairwise_termsim(x) %>%
        treeplot(showCategory = 50, offset.params = list(extend = ratio))
}

remove_cellcycle <- function(
        x,
        regex = c(
            "spindle",
            "mitotic",
            "G2",
            "G1",
            "chromati",
            "meiosis",
            "chromosom",
            "DNA",
            "[cC]ell cycle",
            "organelle",
            "cytokinesis",
            "((nuclear)|(cell)) division",
            "tubule",
            "nucleosome",
            "ATP",
            "NAD",
            "kinetochore",
            "double-strand",
            "recombin",
            "CMG",
            "naphase",
            "hromocenter",
            "ronucleus",
            "idbody"
        )
) {
    x@result <- x@result  %>%
        filter(str_detect(Description, paste0(regex, collapse= "|"), negate = TRUE))
    return(x)
}

filter_gsea <- function(x, regex, FDR = 0.05, width = 500) {
    x@result <- x@result %>%
        filter(p.adjust <= 0.05) %>%
        mutate(p.adjust = as.numeric(p.adjust))
    if (!is.null(regex)) {
        x@result <- x@result %>%
            filter(str_detect(Description, regex))
    }
    x@result[, "Description"] <- to_title(x@result[, "Description"]) %>% str_wrap(width)
    return(x)
}

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
        colour = c(GimmeMyPlot:::palette_discrete()[1], "gray50", GimmeMyPlot:::palette_discrete()[2]),
        title_size = "# Leading genes",
        ...
) {
    p <- filter_gsea(x, regex = regex, FDR = FDR, width = width) %>%
        pairwise_termsim(showCategory = showCategory + 1) %>%
        emapplot(cluster.params = cluster.params, cex.params = cex.params, alpha = alpha, showCategory = showCategory, ...) +
        theme_void() +
        theme(legend.title = element_text(face = "italic", size = 12 * cex))

    theme_enrich0(p, cex, colour, num = FALSE, title_size = title_size)
}
