format_deg <- function(
        x,
        fc_threshold = log2(1.5),
        p_threshold = 0.05
) {
    x %>%
        mutate(
            log10p = -log(padj, 10),
            name = rownames(.) %>% str_remove_all("_\\d*"),
            Expression = case_when(
                log2FoldChange >= fc_threshold & padj <= p_threshold ~ "Up-regulated",
                log2FoldChange <= -fc_threshold & padj <= p_threshold ~ "Down-regulated",
                TRUE ~ "ns"
            )) %>%
        mutate(padj = ifelse(padj == 0, min(padj[padj > 0]), padj)) %>%
        filter(log2FoldChange >= 0.01 | log2FoldChange <= -0.01)
}

extract_top <- function(res, fc_threshold = 1, p_threshold = 0.05, n = 1000, rank = FALSE, var = "pfc", f = desc) {
    temp <- res %>%
        mutate(
            pfc = -log10(padj) * log2FoldChange,
            rank_p = rank((log10(padj))),
            rank_fc = rank(desc(abs(log2FoldChange))),
            rank_pfc = rank(rank_p + rank_fc) / 2
        ) %>%
        dplyr::arrange(f(abs(!!sym(var)))) %>%
        filter(abs(log2FoldChange) >= fc_threshold & padj <= p_threshold) %>%
        head(n)
    if (!rank) {
        temp <- dplyr::select(temp, -c(pfc, contains("rank")))
    }
    return(temp)
}

theme_bulk <- function(p, cex = 1, show_axis = TRUE, colour = "black") {
    axis <- element_text(
        # face = "bold.italic",
        size = 18 * cex
    )
    p <- p +
        theme(
            axis.title = axis,
            axis.text = element_text(size = 15 * cex, color = colour),
            axis.ticks = element_line(color = colour, linewidth = 1), #size = 1.2),
            plot.title = element_text(size = 22 * cex, hjust = 0.5),
            legend.title = element_text(size = 12 * cex),
            legend.text = element_text(colour = colour, size = 10 * cex),
            panel.border = element_rect(colour = colour, fill = NA, size = 1)
        )
    if (show_axis) {
        p <- p + theme(axis.line = element_line(color = colour, linewidth = 1)) #size = 1))
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

log0x_trans <- trans_new("log0x", function(x) x, function(x) x, breaks = breaks_extended(4), format = label_number_auto())
log1x_trans <- trans_new("log1x", function(x) logx_trans(x, 1), function(x) expx_trans(x, 1), breaks = breaks_extended(4), format = label_number_auto())
log2x_trans <- trans_new("log2x", function(x) logx_trans(x, 2), function(x) expx_trans(x, 2), breaks = breaks_extended(4), format = label_number_auto())
log10x_trans <- trans_new("log10x", function(x) logx_trans(x, 10), function(x) expx_trans(x, 10), breaks = breaks_extended(4), format = label_number_auto())

log10nx_trans <- trans_new("log10nx", function(x) -log10(x), function(x) 10^(-x), breaks = breaks_extended(6), format = label_number_auto())

format_labels <- function(x) {
    labels <- scales::label_number_auto()(x)
    x <- as.character(x)
    x[x == "0.0"] <- "0"
    x[x == "1e+00"] <- "1"
    return(x)
}

volcano_plot <- function(
        res,
        top_genes = NULL,
        title = "",
        legend = "right",
        cex = 1.5,
        fc_threshold = 0.5,
        p_threshold = 0.05,
        force = 10,
        fc_log = 2,
        width = 30,
        ylab = NULL,
        ds = NULL,
        ps = NULL,
        ...) {
    if (fc_log == 1) {
        log_func <- log1x_trans
    } else if (fc_log == 2) {
        log_func <- log2x_trans
    } else if (fc_log == 10) {
        log_func <- log10x_trans
    } else {
        log_func <- log0x_trans
    }
    # seq(3, 10) %>% paste0("2^", .) %>% sapply(function(x) parse(text = x) %>% eval())
    # ds <- c(1.5, 3, 6, 20, 60)
    if(is.null(ds)) {
        if (is.null(ylab))
            ds <- c(min(res$log2FoldChange, na.rm = TRUE), max(res$log2FoldChange, na.rm = TRUE)) %>% pretty(5) %>% expx_trans(base = fc_log)
        # if (fc_log > 0 && min(res$log2FoldChange, na.rm = TRUE) <= 1) {
        # ds <- c(0.5, 1, 1.25, 1.5, 2, 3, 6, 9, 12)
        # ds <- c(1.5, 2.5, 5, 25, 150)
        else
            ds <- ylab
        # }
    }
    if(is.null(ps))
        ps <- min(res$padj[res$padj!=0], na.rm = TRUE) %>% log10() %>% `-`(1) %>% seq(-1, .) %>% pretty(3) %>% paste0("1e", .) %>% as.numeric() %>% c(1)
    # ps <- c(0, 0.5, 0.2, 0.1)
    # ps <- c(1, 1e-15, 1e-30)
    # ds <- c(-1.2, -0.8, -0.5, 0.5, 0.8, 1.2)
    res <- select(res, log2FoldChange, padj, name, Expression)
    tmp <- data.frame(t(c(Inf, -1, NA))) %>%
        mutate_all(as.numeric) %>%
        data.frame(c("Up-regulated", "Down-regulated")) %>%
        set_colnames(colnames(res))
    res0 <- rbind(tmp, res)
    p <- ggplot(res0, aes(expx_trans(log2FoldChange, fc_log), padj)) +
        geom_point(aes(fill = Expression), size = 3 * cex, alpha = 0.5, shape = 21, stroke = NA) +
        geom_vline(xintercept = expx_trans(fc_threshold, fc_log) * c(-1, 1), colour = "gray50", lty = 2, lwd = 1.2) +
        geom_hline(yintercept = p_threshold, colour = "gray50", lty = 2, lwd = 1.2) +
        xlab("Fold change") +
        ylab("False Discovery Rate") +
        scale_fill_manual(values = c("#A6CEE3", "gray80", "#FB9A99")) +
        scale_color_manual(values = c("dodgerblue3", "firebrick3")) +
        guides(colour = guide_legend(override.aes = list(size = 2))) +
        scale_x_continuous(trans = log_func, breaks = c(rev(ds) * -1, ds) %>% round(1), labels =  function(x) format_labels(x)) +
        scale_y_continuous(breaks = ps, trans = log10nx_trans, labels = function(x) format_labels(x)) +
        ggtitle(label = title) +
        theme_classic() %>%
        theme_bulk(cex) +
        theme(legend.position = legend)
    if (!is.null(top_genes)) {
        top_genes <- select(top_genes, log2FoldChange, padj, name, Expression) %>%
            rbind(tmp)
        p + geom_text_repel(
            data = top_genes,
            mapping = aes(
                expx_trans(log2FoldChange, fc_log),
                padj,
                label = str_wrap(name, width),
                colour = Expression
            ),
            size = cex * 6,
            force = force,
            # fontface = "italic",
            segment.color = "grey50",
            # colour = "black",
            ...
        )
    } else {
        p
    }
}


print_dgea0 <- function(x, base = 2, ...) {
    if (base == 2) {
        func <- function(x) 2^(x)
    } else {
        func <- function(x) exp(x)
    }
    extract_top(x, n = 10000, rank = TRUE, fc_threshold = 0, p_threshold = 1, ...) %>%
        filter(Expression != "ns") %>%
        mutate(
            FC = ifelse(
                log2FoldChange > 0,
                round(func(log2FoldChange), 2),
                -round(func(abs(log2FoldChange)), 2)
            )
        ) %>%
        select(name, contains("ensembl_ids"), baseMean, FC, log2FoldChange, contains("lfcSE"), pvalue, padj, log10p, pfc, rank_fc, rank_p, rank_pfc, Expression)
}

print_dgea <- function(x, name = FALSE, ensembl = NULL, base = 2, ...) {
    res <- print_dgea0(x, base = base, ...)
    if (name) {
        res <- mutate(res, full_name = get_ncbi_name(res, ensembl = ensembl))
    }
    res %>%
        rename(alias = name) %>%
        mutate(FDR = format(padj, digits = 2, scientific = TRUE)) %>%
        select(alias, contains("full_name"), FC, FDR)
}

get_ncbi_description <- function(x, metadata_genes) {
    gene_subset <- filter(metadata_genes, ensembl_gene_id == x)
    gene_entrez <- pull(gene_subset, "entrezgene_id") %>%
        .[1]
    if(!is.na(gene_entrez))
        paste0("https://www.ncbi.nlm.nih.gov/gene/", gene_entrez) %>%
        read_html() %>%
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
        GimmeMyPlot:::to_title() %>%
        paste0(pull(gene_subset, "gene_name"), ": ", .)
    else
        paste0(pull(gene_subset, "gene_name"), ": ", NA)
}

plot_heatmap <- function(x, transpose = FALSE, ...) {
    df_genes <- assay(count_normalized)[pull(x, 1), ] %>%
        set_colnames(pull(metadata_samples, 1)) %>%
        set_rownames(pull(x, gene_name))
    sampleDists <- dist(t(df_genes))
    geneDists <- dist(df_genes)
    plot.new()
    pheatmap(
        as.matrix(df_genes),
        clustering_distance_rows = geneDists,
        clustering_distance_cols = sampleDists,
        col = colorRampPalette(brewer.pal(9, "Reds"))(255),
        ...
    )
}
