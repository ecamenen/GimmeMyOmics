#' Displays Enrichment
#'
#' Formats and filters enrichment analysis results from enrichment analysis libraries:
#' `enrichR` and `clusterProfiler`.
#'
#' @param x An enrichment results.
#' @param method A character string specifying the enrichment analysis method. Options are:
#'   - `"enrichr"`: Enrichr results (**requires `path2gene`**, from `enrichR` package).
#'   - `"gsea"`: Gene Set Enrichment Analysis (GSEA) results (from `clusterProfiler` package).
#'   - `"ora"`: Over-enrichment Analysis (GSEA) results (from `clusterProfiler` package).
#' @param regex A character vector of regular expressions to filter terms by their descriptions.
#' @param pval A numeric threshold for filtering results based on FDR-adjusted p-values.
#'
#' @details
#' This function standardizes different enrichment result formats, applying transformations such as:
#' - Computing the percentage of differentially expressed genes (DEGs) in each pathway.
#' - Formatting gene lists for better readability.
#' - Filtering results based on FDR and optional regex pattern matching.
#'
#' @return A formatted `tibble` containing processed enrichment results,
#'   with columns:
#'   - `Description`: The pathway or term name.
#'   - `FDR`: The adjusted p-value.
#'   - `Nb DEG` / `Nb genes`: Number of differentially expressed genes and total genes in the pathway.
#'   - `DEG/Genes`: Percentage of DEGs in the pathway.
#'   - `Genes`: List of genes involved.
#'   - Additional columns depending on the enrichment method (e.g., `NES`, `ID`).
#'
#' @examples
#' # Example 1: Enrichr results
#' path_name <- c("Neutrophil degranulation", "Macrophage migration")
#' pval <- c(2.7e-02, 2.89e-03)
#' genes <- c("ITGB2/ANXA3/STXBP2/SPI1/ITGAM/CD177", "MAPK3/AKIRIN1/CX3CR1/CNN2/LGALS3/B4GALT1/C3AR1")
#' ids <- c("GO:0043312", "GO:1905517")
#' 
#' enrichr_results <- data.frame(
#'   Term = path_name,
#'   Adjusted.P.value = pval,
#'   Overlap = c("6/12", "7/17"),
#'   ID = ids,
#'   Genes = gsub("/", ";", genes)
#' )
#' print_enrich(enrichr_results, method = "enrichr")
#'
#' # Example 2: GSEA results (from clusterProfiler)
#' gsea_results <- data.frame(
#'   Description = path_name,
#'   p.adjust = pval,
#'   core_enrichment = genes,
#'   setSize = c(12, 17),
#'   NES = c(1.75, 2.01),
#'   ID = ids
#' )
#' print_enrich(gsea_results, method = "gsea")
#'
#' # Example 3: Over-representation enrichment results (from clusterProfiler)
#' ora_results <- data.frame(
#'   Description = path_name,
#'   p.adjust = pval,
#'   GeneRatio = c("6/12", "7/17"),
#'   BgRatio = c("60/120", "70/140"),
#'   geneID = genes,
#'   ID = ids
#' )
#' print_enrich(ora_results)
#'
#' @export
print_enrich <- function(x, method = "ora", regex = NULL, pval = 0.05) {
    x <- as_tibble(x)

    # if (method == "enrichr" && is.null(path2gene)) {
    #     stop("`path2gene` parameter must not be empty for Enrichr.")
    # }
    # 
    # if (!is.null(path2gene)) {
    #     x <- x %>% mutate(
    #         x,
    #         bg = list.mapv(
    #             x[[1]] %>% as.data.frame() %>% pull(1),
    #             f(i) ~
    #                 path2gene[pull(path2gene, 1) %in% i, ] %>%
    #                 pull(2) %>%
    #                 length()
    #         )
    #     )
    # }

    x <- switch(
        method,
        "enrichr" = x %>%
          mutate(
            Description = Term,
            FDR = Adjusted.P.value,
            Count = as.numeric(str_split_fixed(Overlap, "/", 2)[, 1]),
            bg = as.numeric(str_split_fixed(Overlap, "/", 2)[, 2])
          ),
        "gsea" = x %>%
            mutate(
                FDR = p.adjust,
                Count = sapply(stri_split_fixed(core_enrichment, "/"), function(i) length(unique(i))),
                bg = setSize,
                Genes = sapply(stri_split_fixed(core_enrichment, "/"), function(i) paste(unique(i), collapse = ";"))
            ) %>%
            rename(pval = "p.adjust"),
        x %>%
          mutate(
            x,
            FDR = p.adjust,
            Count = str_split_fixed(GeneRatio, "/", 2) %>% .[, 1] %>% as.numeric(),
            bg = str_split_fixed(BgRatio, "/", 2) %>% .[, 1] %>% as.numeric(),
            Genes = str_replace_all(geneID, "/", ";")
        ) %>%
          select(-geneID)
    )

    if (!is.null(regex)) {
        x <- x %>% filter(str_detect(Description, paste(regex, collapse = "|")))
    }

    x <- x %>% mutate(Description = to_title(Description) %>% str_remove(" - .*"))

    res <- x %>%
        filter(FDR <= pval) %>%
        rename(`Nb DEG` = Count, `Nb genes` = bg) %>%
        mutate(
            FDR = format(FDR, digits = 3, scientific = TRUE),
            `DEG/Genes` = round(`Nb DEG` / `Nb genes` * 100, 1)
        )%>%
        select(Description, FDR, `Nb DEG`, `Nb genes`, `DEG/Genes`, Genes, contains(c("NES", "ID")))

    if (method == "gsea") {
        res %>%
            rename_with(~ str_replace_all(., "DEG", "Enriched"), contains("DEG")) %>%
            relocate("NES", .after = "FDR") %>%
            relocate("ID", .after = "Description")
    } else {
      res
    }
}

theme_enrich <- function(
        p,
        cex = 1,
        colour_gradient = brewer.pal(11, "Spectral") %>% rev(),
        colour_text = "black",
        title = NULL,
        label_x = "Gene ratio",
        title_size = "# DEG") {
    p <- p + theme_minimal() %>%
        theme_bulk(cex * 1.5) +
        theme(
            axis.ticks.y = element_blank(),
            axis.text.y = element_text(
                size = 20 * cex,
                colour = colour_text
            ),
            axis.text.x = element_text(
                size = 20 * cex,
                angle = 45,
                hjust = 1,
                vjust = 1,
                colour = "black"
            ),
            axis.title = element_text(
                face = "bold.italic",
                size = 30 * cex
            )
        ) +
        labs(x = label_x, y = NULL) +
        scale_x_continuous(breaks = pretty_breaks(n = 3))
    if (label_x == "Gene ratio")
        p <- p + scale_x_continuous(labels = label_percent(1), breaks = pretty_breaks(n = 3))
    return(theme_enrich0(p, cex, colour_gradient, title, "FDR", title_size = title_size))
}

theme_enrich0 <- function(
        p,
        cex = 1,
        colour = c(palette_discrete()[1], "gray50", palette_discrete()[2]),
        title = NULL,
        label_colour = "FDR",
        range = c(5, 12) * cex,
        trans = FALSE,
        power = 2,
        title_size = "# Leading genes",
        num = TRUE) {
    if (label_colour == "FDR") {
        label_func <- label_pvalue()
    } else {
        label_func <- label_number_auto()
        if (!isFALSE(trans))
            label_func <- function(x) expx_trans(x, base = power) %>% round(1)
    }
    p <- p +
        scale_fill_gradientn(
            labels = function(x) label_func(x),
            breaks = breaks_pretty(n = 4),
            name = label_colour,
            colours = colour
        ) +
        scale_size_continuous(
            range = range,
            breaks = function(x) unique(round(pretty(x, n = 4))),
            labels = label_number(accuracy = 1),
            name = title_size
        ) +
        guides(
            size = guide_legend(
                order = 1,
                override.aes = list(fill = "black", color = "black")
            ),
            fill = guide_colorbar(
                order = 2,
                frame.colour = "black",
                frame.linewidth = 0.75,
                ticks.colour = "black",
                ticks.linewidth = 0.75
            ),
            color = guide_colorbar(
                barheight = 10 * cex,
                order = 2,
                frame.colour = "black",
                frame.linewidth = 0.75 * cex,
                ticks.colour = "black",
                ticks.linewidth = 0.75 * cex
            )
        ) +
        labs(title = title) +
        theme(
            plot.title = element_text(
                size = 30 * cex,
                face = "bold"
            ),
            legend.title = element_text(
                face = "bold.italic",
                size = 20 * cex
            ),
            legend.text = element_text(size = 20 * cex)
        )

    if (num) {
        p + scale_color_gradientn(
            labels = function(x) label_func(x),
            breaks = breaks_pretty(n = 4),
            # breaks = breaks_width(width = 5, offset = 0),
            name = label_colour,
            # limits = c(0, NA),
            colours = colour
        )
    } else {
        p
    }
}

#' Plot Enrichment Analysis Results
#'
#' Visualizes the results of enrichment analysis (e.g., GO, GSEA, KEGG, or limma) using dot plots. The function supports customization of the plot appearance, including colors, labels, and titles.
#'
#' @inheritParams print_enrich
#' @inheritParams network_enrich
#' @param n Integer, the number of top terms to display.
#' @param cex Numeric, the scaling factor for text size.
#' @param ratio Numeric, the ratio for adjusting plot limits.
#' @param width Integer, the maximum width (in characters) for term labels.
#' @param path2gene Optional, a `data.frame` or `tibble` mapping pathways to genes. Required for `method = "kegg"` if `GeneRatio` is not provided.
#' @param colour Character vector of length 3, specifying the colors for the gradient fill of the points.
#' @param label_x Character, the label for the x-axis. Default is `"generatio"` (gene ratio). For `method = "gsea"`, this can be changed to another column name.
#'
#' @details
#' The function performs the following steps:
#' 1. Filters and processes the input data based on the `method` parameter.
#' 2. Computes the gene ratio for each term.
#' 3. Truncates and formats term labels for readability.
#' 4. Ranks terms by significance and selects the top `n` terms.
#' 5. Generates a dot plot with term labels on the y-axis, gene ratio on the x-axis, and point size proportional to the number of genes.
#'
#' For `method = "gsea"`, the function also extracts the number of leading genes from the `core_enrichment` column.
#'
#' @return A `ggplot` object representing the enrichment plot.
#'
#' @examples
#' # Example 1: Enrichr results (Requires a named list as path2gene)
#' path_name <- c("Neutrophil degranulation", "Macrophage migration")
#' pval <- c(2.7e-02, 2.89e-03)
#' genes <- c("ITGB2/ANXA3/STXBP2/SPI1/ITGAM/CD177", "MAPK3/AKIRIN1/CX3CR1/CNN2/LGALS3/B4GALT1/C3AR1")
#' ids <- c("GO:0043312", "GO:1905517")
#'
#' enrichr_results <- data.frame(
#'   Term = path_name,
#'   Adjusted.P.value = pval,
#'   Overlap = c("6/12", "7/17"),
#'   ID = ids,
#'   Genes = gsub("/", ";", genes)
#' )
#' plot_enrich(enrichr_results, method = "enrichr")
#'
#' # Example 2: GSEA results (from clusterProfiler)
#' gsea_results <- data.frame(
#'   Description = path_name,
#'   p.adjust = pval,
#'   core_enrichment = genes,
#'   setSize = c(12, 17),
#'   NES = c(1.75, 2.01),
#'   ID = ids
#' )
#' plot_enrich(gsea_results, method = "gsea")
#'
#' # Example 3: Over-representation enrichment results (from clusterProfiler)
#' ora_results <- data.frame(
#'   Description = path_name,
#'   p.adjust = pval,
#'   GeneRatio = c("6/12", "7/17"),
#'   BgRatio = c("60/120", "70/140"),
#'   geneID = genes,
#'   ID = ids
#' )
#' plot_enrich(ora_results)
#' 
#' # Example GO enrichment results
#' go_results <- data.frame(
#'   Description = c("immune response", "cell cycle", "DNA repair"),
#'   p.adjust = c(0.001, 0.01, 0.05),
#'   GeneRatio = c("10/100", "15/150", "20/200"),
#'   BgRatio = c("100/1000", "150/1500", "200/2000")
#' )
#'
#' @export
plot_enrich <- function(
        x,
        n = 20,
        title = NULL,
        method = "ora",
        cex = 0.65,
        ratio = 5,
        width = 50,
        path2gene = NULL,
        colour = c(palette_discrete()[1], "grey80", palette_discrete()[2]),
        regex = NULL,
        label_x = "generatio") {
    func <- function(x) {
        str_remove_all(x, "Genes ((down)|(up))-regulated ((in ?)|(with))") %>%
            str_remove_all("comparison of ") %>%
            str_remove_all("Genes ((posi)|(nega))tively correlated with ") %>%
            str_remove_all("\\[GeneID=\\d*\\]") %>%
            str_remove_all("([uU]ntreated )?peripheral blood mono((nuclear)|(cytes))( cells?)?( \\(PBMC\\))?( from)? ") %>%
            str_remove_all("the ") %>%
            str_remove(" - .*")
    }
    if (method == "gsea") {
        if (!is.null(regex)) {
            x <- filter_gsea(x, regex)
        }
        df <- mutate(
            x,
            Adjusted.P.value = p.adjust,
            Term = Description %>%
                func() %>%
                to_title(),
            Count = str_split(core_enrichment, "/") %>% sapply(length),
            Overlap = core_enrichment
        )
        title_size <- "# Leading genes"
        if (label_x == "generatio") {
            df <- mutate(df, generatio = Count / setSize)
            label_x <- "Gene ratio"
        } else {
            df <- mutate(df, generatio = !!sym(label_x))
        }
        df <- arrange(df, desc(generatio))
        x_var <- "Adjusted.P.value"
    } else if (method == "ora") {
        title_size <- "# DEG"
        label_x <- "Gene ratio"
        if (!is.null(regex)) {
          x <- filter_gsea(x, regex)
        }
        df <- mutate(
            x,
            Adjusted.P.value = p.adjust,
            Term = Description %>%
                func() %>%
                to_title(),
            Count = str_split(genes, "/") %>% sapply(length)
        )
        if (!is.null(path2gene)) {
            n_paths <- list.mapv(
                pull(x, 1),
                f(i) ~
                    path2gene[pull(path2gene, 1) %in% i, ] %>%
                    pull(2) %>%
                    length()
            )
            df <- mutate(df, bg = n_paths, generatio = Count / n_paths)
        } else {
            df <- mutate(df, bg = n, generatio = str_split(GeneRatio, "/") %>%
                             sapply(function(i) as.numeric(i[1]) / as.numeric(i[2]))
            )
        }
        x_var <- "p.adjust"
    } else {
        title_size <- "# DEG"
        label_x <- "Gene ratio"
        if (!is.null(regex)) {
            x <- filter(x, str_detect(Term, paste(regex, collapse = "|")))
        }
        df <- mutate(
            x,
            Count = str_remove_all(Overlap, "\\/.*") %>% as.numeric(),
            generatio = {
                str_split(Overlap, "/") %>%
                    sapply(function(i) as.numeric(i[1]) / as.numeric(i[2]))
            },
            bg = str_remove_all(Overlap, ".*\\/") %>% as.numeric(),
            Term = Term %>%
              func() %>%
              to_title()
        )
        x_var <- "Adjusted.P.value"
    }
    df0 <- filter(df, Adjusted.P.value <= 0.05) %>%
        filter(!is.na(Term))
    if (nrow(df0) < n)
        df0 <- df
    df0 <- arrange(df0, abs(!!sym(x_var)))
    y <- "generatio"
    # if (method %in% c("gsea", "kegg")) {
    #   df0 <- arrange(df0, Adjusted.P.value)
    #   y <- "Adjusted.P.value"
    # } else {
    #   df0  <- arrange(df0, desc(Combined.Score))
    #   y <- "Combined.Score"
    # }
    df <- head(df0, n) %>%
        arrange(!!sym(x_var)) %>%
        mutate(
            label = {
                str_remove_all(Term, "\\(.*\\)") %>%
                    str_remove_all("((ORPHA)|(WP)|(HSA)|(R-)|(CL:)).*") %>%
                    str_pretty(width) %>%
                    str_trim() %>%
                    to_title()
            },
            rank = rev(row_number(!!sym(y)))
        )
    # if (method %in% c("gsea", "kegg")) {
    #   df <- mutate(df, rank = rev(row_number(!!sym(y))))
    # } else {
    df <- mutate(df, rank = row_number(!!sym(y)))
    # }
    # print(as_tibble(df) %>% select(1, 2, 4))
    if (method == "gsea") {
        colour_path <- "black"
    } else {
        colour_path <- ifelse(df$Adjusted.P.value <= 0.05, palette_discrete()[1], "gray50")
    }
    p <- ggplot(df, aes(generatio, rank)) +
        geom_point(
            aes(fill = Adjusted.P.value, size = Count),
            colour = "black",
            pch = 21,
            # stroke = NA
        ) +
        scale_y_continuous(breaks = df$rank, labels = df$label)

    theme_enrich(
        p,
        cex,
        colour_gradient = colour,
        colour_text = colour_path,
        title = title,
        label_x = label_x,
        title_size = title_size
    )
    # expand_limits(y = max(df$generatio) + max(df$generatio) / ratio)
}

#' Pathway keyword dictionary for immune-related terms
#' 
#' Creates a categorized list of regular expression patterns for identifying
#' immune-related pathways in enrichment analyses. The patterns are organized
#' hierarchically from specific to broad immune categories.
#'
#' @return A nested list containing regular expression patterns for immune-related
#' pathway identification. The list contains the following categories:
#' \describe{
#'   \item{cell}{Immune cell types (neutrophils, macrophages, T cells, B cells, etc.)}
#'   \item{cytokine}{Cytokines and lipid mediators (interleukins, interferons, prostaglandins, etc.)}
#'   \item{cytokine_full}{Extended cytokine patterns including TNF and NF-κB}
#'   \item{immunity}{Core immune system terms (cells + cytokines + basic immune processes)}
#'   \item{immunity_additional}{Extended immunity terms including complement system and TLRs}
#'   \item{immunity_full}{Comprehensive immune system terms including inflammation and hematopoiesis}
#' }
#'
#' @examples
#' \dontrun{
#' # Get all immune-related keywords
#' keywords <- pathway_keywords()
#' 
#' # Search for macrophage-related pathways
#' macrophage_pathways <- str_detect(
#'   pathway_descriptions,
#'   keywords$cell[2]
#' )
#' 
#' # Create a comprehensive immune filter
#' immune_filter <- paste(unlist(keywords$immunity_full), collapse = "|")
#' immune_pathways <- str_subset(pathway_descriptions, immune_filter)
#' }
#'
#' @export
pathway_keywords <- function() {
    l <- list()
    l[["cell"]] <- c("eutrophil", "(acrophage)|(onocyte)", "endritic cell", "((natural killer)|(NK)) cell", "(T [- ]? cell)|(T-helper)|(CD[48][- ])", "B[- ]?cell", "NETosis", "Th\\d{1,2} cell")
    l[["cytokine"]] <- c("(nterleukins?)|(IL-?\\d{1,2})", "(nterferon)|(IFN[ABG])", "(rostaglandin)|([Aa]rachidonic)|(icosa)|([Ll]eukotriene)|([Dd]ocosahexaenoic)|([Ii]cosapentaenoic)|([Ll]ipoxin)|(esolvin)")
    l[["cytokine_full"]] <- c(l[["cytokine"]], "(tumor necrosis factor)|(TNF)|(NF-k)")
    l[["immunity"]]  <- c(l[["cell"]], l[["cytokine"]], "STAT[ 35]", "AGE", "[Ll]upus", "(steo[cb]last)|([Bb]one)|(keletal)|(ossification)", "[Aa]rthrit", "[Gg]lucocorticoid")
    l[["immunity_additional"]] <- c(l[["immunity"]], l[["cytokine_full"]], "(omplement)|([^ ]C2 )", "(oll-like)|(TLR )", "mTORC1", "(Fc gamma)|(FCG)", "etalloproteinas", "[Aa]cute") %>% unique()
    l[["immunity_full"]] <- c(l[["immunity_additional"]],  "(PUMA)|(TP53)|( p53)", "([Ii]nflamm)|([Ii]mmun)", "hemokine", "mhc",  "phago((cytosis)|(some))", "leukocyte", "myeloid", "cytokine[^sis]", "granulocyte", "[Ll]ympho", "[Hh]emopo")
    l[["cell_cycle"]] <- c(
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
    return(l)
}

format_path <- function(x, width = 20) {
    str_remove_all(x, "\\(.*\\)") %>%
        str_remove_all("((ORPHA)|(WP)|(HSA)|(R-)|(CL:)).*") %>%
        # str_remove_all("(Signaling)|(By)") %>%
        # str_replace_all("Interleukin", "IL") %>%
        # str_replace(" [aA]nd ", "\\/")
        str_pretty(width) %>%
        str_trim() %>%
        # sort() %>%
        to_title()
}

list_table <- function(x) {
    ids <- unlist(x) %>% unique()
    list.map(
        ids,
        f(i) ~ list.if(x, i %in% .) %>% as.numeric()
    ) %>%
        list.rbind() %>%
        set_colnames(names(x)) %>%
        as.data.frame() %>%
        select(colnames(.))
}

list_common <- function(x) {
    ids <- unlist(x) %>% unique() %>% sort()
    res <- list()
    for (i in ids) {
        tmp <- list.which(x, i %in% .)
        if (length(tmp) > 1) {
            n <- paste(names(x)[tmp] %>% sort(), collapse = ":")
            if(n %in% names(res)) {
                res[[n]] <- c(res[[n]], i)
            } else {
                res[[n]] <- i
            }
        }
    }
    res <- res[order(names(res))]
    len0 <- names(res) %>% str_extract_all(":") %>% sapply(length) %>% set_names(seq(.)) %>% sort()
    len <- names(len0) %>% as.numeric()
    return(res[len])
}

list_count <- function(x) {
    list_common(x) %>%
        list.class(
            f(x, y, z) ~ z %>%
                list.map(
                    f(x) ~str_split(x, ":") %>%
                        pluck(1) %>%
                        length()
                )
        )
}

heatmap_enrich <- function(x, cex = 1, width_text = 20, power = 2) {
    x %>%
        mutate(ID = rownames(.)) %>%
        gather("key", "value", -ID) %>%
        mutate(
            key = str_wrap(key, width_text),
            ID = str_wrap(ID, width_text)
        ) %>%
        left_join(gene_path1, by = "ID") %>%
        mutate(
            value2 = ifelse(value == 0, NA, value2),
            key = factor(key, levels = colnames(x) %>% str_wrap(width_text))
        ) %>%
        ggplot(aes(ID, key, fill = value2)) +
        geom_tile() +
        scale_fill_gradientn(
            colours = brewer.pal(11, "Spectral") %>% rev(),
            na.value = "white",
            name = "Fold Change",
            breaks = pretty_breaks(5),
            labels = function(x) expx_trans(x, base = power) %>% round(1)
            # name = "-log10(FDR) \n* log2(FC)"
        ) +
        theme_classic() +
        theme(
            axis.title = element_blank(),
            axis.text.x = element_text(size = 10 * cex, color = "grey40", angle = 90, hjust = 1, vjust = 0.5),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            axis.text.y = element_text(size = 10 * cex, color = "grey40"),
            legend.title = element_text(face = "bold.italic", size = 12 * cex),
            legend.text = element_text(size = 10 * cex, color = "grey40")
        )
}
