#' Format Differential Expression Analysis (DEA) Results
#'
#' Formats the results of a differential expression analysis (DEA) by adding relevant columns, filtering, and categorizing genes based on fold change and adjusted p-value thresholds.
#'
#' @param x A `data.frame` or `tibble` containing DEA results. Expected columns include `log2FoldChange`, `padj`, and row names representing Ensembl gene identifiers.
#' @param metadata_genes A `data.frame` containing gene metadata, including the columns `ensembl_gene_id` and `gene_name`.
#' @param fc_threshold Numeric, the absolute log2 fold change threshold. Default is `log2(2)`.
#' @param p_threshold Numeric, the adjusted p-value threshold. Default is `0.05`.
#'
#' @details
#' The function performs the following steps:
#' 1. Adds a `log10p` column, which is the -log10 of the adjusted p-value (`padj`).
#' 2. Merges gene metadata to include gene names.
#' 3. Categorizes genes into "Up-regulated," "Down-regulated," or "ns" (not significant) based on the provided fold change and p-value thresholds.
#' 4. Filters out genes with very small fold changes (absolute `log2FoldChange` < 0.01).
#'
#' @return A `tibble` with the following additional columns:
#' \describe{
#'   \item{log10p}{Numeric, the -log10 of the adjusted p-value.}
#'   \item{gene_name}{Character, the official gene name from HGNC.}
#'   \item{Expression}{Factor, the expression category ("Up-regulated," "Down-regulated," or "ns").}
#' }
#'
#' @examples
#' # Example DEA results
#' ensembl_gene_id <- c("ENSG00000141510", "ENSG00000012048", "ENSG00000146648",
#'  "ENSG00000157764", "ENSG00000133703")
#' dea_results <- data.frame(
#'   log2FoldChange = c(2.0, -1.8, 0.5, -0.2, 3.0),
#'   padj = c(0.001, 0.01, 0.1, 0.06, 0.0001)
#' )
#' rownames(dea_results) <- ensembl_gene_id
#'
#' # Example gene metadata
#' metadata_genes <- data.frame(
#'   ensembl_gene_id = ensembl_gene_id,
#'   gene_name = c("TP53", "BRCA1", "EGFR", "BRAF", "KRAS")
#' )
#'
#' # Format DEA results
#' format_dea(dea_results, metadata_genes)
#'
#' # Custom thresholds
#' format_dea(dea_results, metadata_genes, fc_threshold = log2(1.5), p_threshold = 0.01)
#' @export
format_dea <- function(
    x,
    metadata_genes = NULL,
    fc_threshold = 1,
    p_threshold = 0.05
) {
    # required_cols <- c("ensembl_gene_id", "gene_name")
    # if (!all(required_cols %in% colnames(metadata_genes))) {
    #     stop("metadata_genes must contain the columns: ", paste(required_cols, collapse = ", "))
    # }
    tmp <- x %>%
        as.data.frame() %>%
        rownames_to_column("ensembl_gene_id")

    if (!is.null(metadata_genes)) {
        tmp <- left_join(
            tmp,
            metadata_genes,
            by = "ensembl_gene_id"
        )
    }
        mutate(
            tmp,
            log10p = -log10(padj),
            gene_name = str_remove_all(gene_name, "_\\d+"),
            Expression = case_when(
                log2FoldChange >= fc_threshold & padj <= p_threshold ~ "Up-regulated",
                log2FoldChange <= -fc_threshold & padj <= p_threshold ~ "Down-regulated",
                TRUE ~ "ns"
            ),
            padj = ifelse(padj == 0, min(padj[padj > 0], na.rm = TRUE), padj)
        ) %>%
        filter(abs(log2FoldChange) >= 0.01) %>%
        relocate("gene_name") %>%
        as_tibble()
}

#' Filter Top Differentially Expressed Genes
#'
#' Filters the top differentially expressed genes from a DEA result set based on fold change, adjusted p-value, and a combined ranking metric.
#'
#' @param x A `data.frame` or `tibble` containing DEA results. Must include `log2FoldChange` and `padj` columns.
#' @param fc_threshold Numeric, minimum absolute log2 fold change required for inclusion.
#' @param p_threshold Numeric, maximum adjusted p-value allowed.
#' @param n Integer, number of top genes to return. Default is `1000`.
#' @param rank Logical, whether to return ranking columns (`rank_p`, `rank_fc`, `rank_pfc`).
#' @param var Character, ranking variable.
#' @param f Function, sorting function (e.g., `desc` for descending).
#'
#' @details
#' The function applies the following steps:
#' - Computes `pfc`, a ranking metric defined as `-log10(padj) * log2FoldChange`.
#' - Calculates separate rankings for significance (`rank_p` from `log10(padj)`) and effect size (`rank_fc` from `abs(log2FoldChange)`).
#' - Computes an average rank (`rank_pfc`) combining both ranks.
#' - Filters genes based on the given fold change and p-value thresholds.
#' - Selects the top `n` genes sorted by the specified ranking variable.
#'
#' @return A `tibble` containing selected genes with additional computed ranking columns if `rank = TRUE`.
#'
#' @examples
#' # Example DEA results
#' dea_results <- data.frame(
#'   gene_name = c("TP53", "BRCA1", "EGFR", "BRAF", "KRAS"),
#'   log2FoldChange = c(2.0, -1.8, 0.5, -0.2, 3.0),
#'   padj = c(0.001, 0.01, 0.1, 0.06, 0.0001)
#' )
#'
#' # Select top genes
#' top_genes(dea_results)
#'
#' # Custom thresholds and return rankings
#' top_genes(dea_results, fc_threshold = log2(1.5), p_threshold = 0.01, 
#' rank = TRUE, var = "log2FoldChange")
#'
#' @export
top_genes <- function(
    x,
    fc_threshold = 1,
    p_threshold = 0.05,
    n = Inf,
    rank = FALSE,
    var = "pfc",
    f = desc
) {
    res <- x %>%
        mutate(
            pfc = -log10(padj) * log2FoldChange,
            rank_p = dense_rank(log10(padj)),
            rank_fc = dense_rank(desc(abs(log2FoldChange))),
            rank_pfc = dense_rank((rank_p + rank_fc) / 2)
        ) %>%
        arrange(f(abs(.data[[var]]))) %>%
        filter(abs(log2FoldChange) >= fc_threshold, padj <= p_threshold) %>%
        slice_head(n = n)

    if (!rank) {
        res <- select(res, -c(pfc, starts_with("rank")))
    }

    return(res)
}


theme_bulk <- function(p, cex = 1, show_axis = TRUE, colour = "black") {
    axis <- element_text(size = 18 * cex)
    p <- p +
        theme(
            axis.title = axis,
            axis.text = element_text(size = 15 * cex, color = colour),
            axis.ticks = element_line(color = colour, linewidth = 1),
            plot.title = element_text(size = 22 * cex, hjust = 0.5),
            legend.title = element_text(size = 12 * cex),
            legend.text = element_text(colour = colour, size = 10 * cex),
            panel.border = element_rect(colour = colour, fill = NA, size = 1)
        )
    if (show_axis) {
        p <- p + theme(axis.line = element_line(color = colour, linewidth = 1))
    }
    p
}

expx_trans <- function(x, base = 2) {
    if (base == 1) {
        exp(abs(x)) * sign(x)
    } else if (base > 1) {
        base^(abs(x)) * sign(x)
    } else {
        return(x)
    }
}

logx_trans <- function(x, base = 2) {
    if (base == 0) {
        return(x)
    }
    if (base == 1) {
        base <- exp(1)
    }
    log(abs(x), base) * sign(x)
}

format_labels <- function(x) {
    labels <- label_number_auto()(x)
    x <- as.character(x)
    x[x == "0.0"] <- "0"
    x[x == "1e+00"] <- "1"
    return(x)
}

#' Volcano Plot for Differential Expression Analysis (DEA) Results
#'
#' Generates a volcano plot to visualize differentially expressed genes based on log fold change and adjusted p-value.
#'
#' @param x A `data.frame` or `tibble` containing DEA results. Must include `log2FoldChange`, `padj`, `gene_name`, and `Expression` columns.
#' @param top_genes Optional. A subset of `x` containing the most significant genes to be labeled on the plot.
#' @param title Optional. A character string specifying the title of the plot.
#' @param legend Character, position of the legend. Default is `"right"`. Other options include `"top"`, `"bottom"`, `"left"`, or `"none"`.
#' @param cex Numeric, scaling factor for text and points.
#' @param fc_threshold Numeric, absolute log2 fold change threshold to determine significance.
#' @param p_threshold Numeric, adjusted p-value threshold to determine significance.
#' @param force Numeric, force applied for label repulsion in `geom_text_repel()`.
#' @param fc_log Numeric, base for log transformation of fold change axis. If `0`, no transformation is applied.
#' @param width Numeric, width for wrapping gene labels in `geom_text_repel()`.
#' @param breaks_x Optional. Numeric vector of breakpoints for the x-axis (fold change). If `NULL`, it is automatically computed.
#' @param breaks_y Optional. Numeric vector of breakpoints for the y-axis (adjusted p-value). If `NULL`, it is automatically computed.
#' @param fontface Character, font style for gene labels (`"plain"`, `"italic"`, `"bold"`, etc.).
#' @param ... Additional arguments passed to `geom_text_repel()`.
#'
#' @details
#' The volcano plot represents genes based on:
#' - **x-axis (Fold change):** Log-transformed fold change (can be in log base 2 or another specified base).
#' - **y-axis (False Discovery Rate - FDR):** Adjusted p-value in a -log10 scale.
#' - **Color coding:** Up-regulated (red), down-regulated (blue), and non-significant (gray) genes.
#' - **Dashed lines:** Represent fold change (`fc_threshold`) and p-value (`p_threshold`) thresholds.
#'
#' If `top_genes` is provided, labels will be added for those genes.
#'
#' @return A `ggplot2` object representing the volcano plot.
#'
#' @examples
#' # Example DEA results
#' dea_results <- data.frame(
#'   gene_name = c("TP53", "BRCA1", "EGFR", "BRAF", "KRAS"),
#'   log2FoldChange = c(2.0, -1.8, 0.5, -0.2, 3.0),
#'   padj = c(0.001, 0.01, 0.1, 0.06, 0.0001),
#'   Expression = c("Up-regulated", "Down-regulated", "ns", "ns", "Up-regulated")
#' )
#'
#' # Generate volcano plot
#' volcano_plot(dea_results)
#'
#' # Highlight top genes
#' tops <- top_genes(dea_results)
#' volcano_plot(dea_results, top_genes = tops)
#'
#' @export
volcano_plot <- function(
        x,
        top_genes = NULL,
        title = NULL,
        legend = "right",
        cex = 1.5,
        fc_threshold = 0.5,
        p_threshold = 0.05,
        force = 10,
        fc_log = 2,
        width = 30,
        breaks_x = NULL,
        breaks_y = NULL,
        fontface = "italic",
        ...
) {
    log_func <- if (fc_log == 0) identity else function(x) logx_trans(x, fc_log)
    exp_func <- if (fc_log == 0) identity else function(x) expx_trans(x, fc_log)

    if (is.null(breaks_x)) {
        breaks_x <- c(min(x$log2FoldChange, na.rm = TRUE), max(x$log2FoldChange, na.rm = TRUE)) %>%
            pretty(5) %>%
            expx_trans(base = fc_log)
    }

    if (is.null(breaks_y)) {
        breaks_y <- min(x$padj[x$padj!=0], na.rm = TRUE) %>%
            log10() %>%
            `-`(1) %>%
            seq(-1, .) %>%
            pretty(3) %>%
            paste0("1e", .) %>%
            as.numeric() %>%
            c(1)
    }

    x <- select(x, log2FoldChange, padj, gene_name, Expression)
    tmp <- data.frame(t(c(.Machine$double.xmax, 1, NA))) %>%
        mutate_all(as.numeric) %>%
        data.frame(c("Up-regulated", "Down-regulated", "ns")) %>%
        set_colnames(colnames(x))

    p <- rbind(x, tmp) %>%
        ggplot(aes(expx_trans(log2FoldChange, fc_log), y = padj)) +
        geom_point(aes(fill = Expression), size = 3 * cex, alpha = 0.5, shape = 21, stroke = NA) +
        geom_vline(xintercept = expx_trans(fc_threshold, fc_log) * c(-1, 1), color = "gray50", linetype = "dashed", lwd = 1.2) +
        geom_hline(yintercept = p_threshold, color = "gray50", linetype = "dashed", lwd = 1.2) +
        scale_fill_manual(values = c("Up-regulated" = "#FB9A99", "Down-regulated" = "#A6CEE3", "ns" = "gray80")) +
        scale_x_continuous(trans = trans_new("lognx", log_func, exp_func), breaks = c(rev(breaks_x) * -1, breaks_x) %>% round(1), labels = format_labels) +
        scale_y_continuous(trans = trans_new("log10x", function(x) -log10(x), function(x) 10^(-x)), breaks = breaks_y, labels = format_labels) +
        xlab("Fold change") +
        ylab("False Discovery Rate") +
        ggtitle(title) +
        theme_classic() %>%
        theme_bulk(cex) +
        theme(legend.position = legend)

    if (!is.null(top_genes)) {
        p <- p +
            scale_colour_manual(values = c("Up-regulated" = "firebrick3", "Down-regulated" = "dodgerblue3", "ns" = "gray80")) +
            guides(colour = "none") +
            geom_text_repel(
            data = top_genes,
            mapping = aes(
                expx_trans(log2FoldChange, fc_log),
                padj,
                label = str_wrap(gene_name, width),
                colour = Expression
            ),
            size = cex * 6,
            force = force,
            segment.color = "grey50",
            fontface = fontface,
            ...
        )
    }

    return(p)
}

#' Print Differential Expression Analysis (DEA) Results
#'
#' Formats and prints the results of a differential expression analysis (DEA), including fold changes, rankings, and additional metadata. The function supports filtering, ranking, and adding gene names from an Ensembl annotation table.
#'
#' @param x A `data.frame` or `tibble` containing DEA results. Expected columns include:
#'   - `gene_name`: Character, the gene symbol or identifier.
#'   - `log2FoldChange`: Numeric, the log2 fold change.
#'   - `padj`: Numeric, the adjusted p-value (FDR).
#'   - `Expression`: Factor, the expression category ("Up-regulated", "Down-regulated", or "ns").
#' @param base Numeric, the base for fold change calculation. Default is `2` (log2 fold change). Use `exp` for natural log fold change.
#' @param ... Additional arguments passed to `top_genes()` for filtering and ranking.
#'
#' @details
#' This function extracts top genes, applies a fold change transformation, filters out non-significant genes, adds gene names (if requested), and formats the FDR values.
#'
#' @return A formatted `tibble` with the following columns:
#' \describe{
#'   \item{alias}{Character, the gene symbol or identifier.}
#'   \item{full_name}{Character, the full gene name (if `name = TRUE` and `metadata_genes` is provided).}
#'   \item{FC}{Numeric, the fold change (transformed based on `base`).}
#'   \item{FDR}{Character, the formatted adjusted p-value (FDR).}
#' }
#'
#' @examples
#' # Example DEA results
#' ensembl_gene_id <- c("ENSG00000141510", "ENSG00000012048", "ENSG00000146648",
#'  "ENSG00000157764", "ENSG00000133703")
#' dea_results <- data.frame(
#'   log2FoldChange = c(2.0, -1.8, 0.5, -0.2, 3.0),
#'   padj = c(0.001, 0.01, 0.1, 0.06, 0.0001)
#' )
#' rownames(dea_results) <- ensembl_gene_id
#'
#' # Example gene metadata
#' metadata_genes <- data.frame(
#'   ensembl_gene_id = ensembl_gene_id,
#'   gene_name = c("TP53", "BRCA1", "EGFR", "BRAF", "KRAS"),
#'   description = c("tumor protein p53", "BRCA1 DNA repair associated",
#'   "KRAS proto-oncogene, GTPase", "epidermal growth factor receptor", 
#'   "B-Raf proto-oncogene")
#' )
#'
#' # Format DEA results
#' dea_results <- format_dea(dea_results, metadata_genes = metadata_genes)
#'
#' print_dea(dea_results)
#'
#' @export
print_dea <- function(x, base = 2, ...) {
    func <- if (base == 2) function(x) 2^x else exp
    top_genes(x, n = 10000, fc_threshold = 0, p_threshold = 1, ...) %>%
    filter(Expression != "ns") %>%
    mutate(
        FC = ifelse(
            log2FoldChange > 0,
            round(func(log2FoldChange), 2),
            -round(func(abs(log2FoldChange)), 2)
        ),
        FDR = format(padj, digits = 2, scientific = TRUE)
    ) %>%
    rename(
        alias = "gene_name",
        full_name = "description"
    ) %>%
    select(alias, full_name, FC, FDR)
}

#' Fetches Gene Description from NCBI
#'
#' Fetches and processes the gene description from NCBI using either an Ensembl or Entrez gene ID.
#'
#' @param x Character, a gene identifier (either Ensembl or Entrez).
#' @param metadata_genes A `data.frame` or `tibble` containing gene metadata with at least `ensembl_gene_id` and `gene_name` columns.
#' @param database Character, specifying the type of gene identifier provided in `x`. Accepted values are `"ensembl_gene_id"` or `"entrezgene_id"`.
#'
#' @details
#' The function performs the following steps:
#' - Constructs the NCBI gene URL and retrieves the gene summary from the webpage.
#' - Cleans the retrieved text by removing redundant phrases and formatting inconsistencies.
#' - Returns a formatted gene description prefixed by the gene symbol.
#'
#' @return A character string containing the formatted gene description, or `NA` if no valid gene ID is found.
#'
#' @examples
#' metadata_genes <- tibble(
#'   ensembl_gene_id = c("ENSG00000141510", "ENSG00000272398"),
#'   entrezgene_id = c(7157, NA),
#'   gene_name = c("TP53", "Unknown Gene")
#' )
#'
#' ncbi_description("ENSG00000141510", metadata_genes)
#' ncbi_description("7157", metadata_genes, database = "entrezgene_id")
#'
#' @export
ncbi_description <- function(x, metadata_genes, database = "ensembl_gene_id") {
    if (!database %in% c("ensembl_gene_id", "entrezgene_id")) {
        stop('Invalid value for `database`. It must match `metadata_genes` column names. Choose either "ensembl_gene_id" or "entrezgene_id".')
    }

    gene_subset <- filter(metadata_genes, .data[[database]] == x) %>%
        slice(1)
    gene_query <- pull(gene_subset, database)

    if (!is.na(gene_query)) {
        url <- paste0("https://www.ncbi.nlm.nih.gov/gene/?term=", gene_query)
        page <- tryCatch(read_html(url), error = function(e) return(NA))

        if (is.na(page)) {
            return(paste0(pull(gene_subset, "gene_name"), ": Website unavailable"))
        }

        page %>%
            html_element("#summaryDl") %>%
            html_text2() %>%
            str_extract("Summary\n(.*)\n", group = 1) %>%
            str_trim() %>%
            str_remove_all("^(T(his)|(he)) ") %>%
            str_remove_all("^((gene)|(locus)|(protein)) ") %>%
            str_remove_all("^encode[sd] (a (matrix )?protein )?((which binds)|(that is))?") %>%
            str_remove_all("^(is) ") %>%
            str_remove_all("^(a member of )*(the)?") %>%
            str_remove_all("^((by )?this gene,?) ") %>%
            str_remove_all("^(that plays a role in) ") %>%
            str_remove_all("^predicted (to )?") %>%
            str_remove_all("^(enables?) (to )?") %>%
            str_remove_all("^(was) ") %>%
            str_remove_all("^belongs? (to )?") %>%
            str_remove_all(". \\[.*\\]$") %>%
            str_trim() %>%
            to_title() %>%
            paste0(pull(gene_subset, "gene_name"), ": ", .)
    } else {
        return(paste0(pull(gene_subset, "gene_name"), ": Description unavailable"))
    }
}

#' Plot Heatmap of Gene Expression Data
#'
#' Plots a heatmap of gene expression data from normalized counts, with optional transposition and clustering.
#'
#' @param x A `data.frame` or `tibble` containing gene identifiers in the first column.
#' @param metadata_samples A `data.frame` or `tibble` containing metadata for the samples.
#' @param count_normalized A matrix or `SummarizedExperiment` object containing normalized gene expression counts.
#' @param transpose Logical, whether to transpose the heatmap.
#' @param colour A color palette used for the heatmap.
#' @param ... Additional arguments passed to `pheatmap`.
#'
#' @details
#' The function performs the following steps:
#' - Extracts expression data using gene identifiers from `x`.
#' - Sets row names to gene names and column names to sample identifiers.
#' - Computes sample and gene distance matrices for clustering.
#' - Plots a heatmap using the `pheatmap` package.
#'
#' @return A heatmap visualization of gene expression.
#'
#' @examples
#' library(DESeq2)
#' # Create example counts matrix
#' nrows <- 5
#' ncols <- 3
#' counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows, ncols)
#' gene_name <- c("TP53", "BRCA1", "EGFR", "BRAF", "KRAS")
#'
#' # Create SummarizedExperiment object
#' count_normalized <- SummarizedExperiment(
#'     assays = SimpleList(counts = counts)
#' )
#' rownames(count_normalized) <- gene_name
#'
#' # Example DEA results
#' dea_results <- data.frame(
#'     gene_name = gene_name,
#'     log2FoldChange = c(2.0, -1.8, 0.5, -0.2, 3.0),
#'     padj = c(0.001, 0.01, 0.1, 0.06, 0.0001)
#' )
#'
#' # Metadata for samples
#' metadata_samples <- data.frame(sample_id = paste0("Sample", 1:ncols))
#' tops <- top_genes(dea_results)
#' tops$ensembl_gene_id <- tops$gene_name
#'
#' # Select top genes and plot heatmap
#' plot_heatmap(tops, metadata_samples, count_normalized)
#'
#' # Transposed heatmap
#' plot_heatmap(tops, metadata_samples, count_normalized, transpose = TRUE)
#'
#' @export
plot_heatmap <- function(
    x,
    metadata_samples,
    count_normalized,
    transpose = FALSE,
    colour = colorRampPalette(brewer.pal(9, "Reds"))(255),
    ...
) {
    gene_counts <- assay(count_normalized)[pull(x, "ensembl_gene_id"), ] %>%
    set_colnames(pull(metadata_samples, 1)) %>%
    set_rownames(pull(x, "gene_name"))

    if (transpose) {
        gene_counts <- t(gene_counts)
    }

    sampleDists <- dist(t(gene_counts))
    geneDists <- dist(gene_counts)

    plot.new()
    pheatmap(
        as.matrix(gene_counts),
        clustering_distance_rows = geneDists,
        clustering_distance_cols = sampleDists,
        color = colour,
        ...
    )
}
