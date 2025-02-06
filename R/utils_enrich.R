library(stringi)

#' @export
#' Cuts a sentence to a given number of characters and leaves the words as integers
#' @example str_trunc1("Hi there, I'm a sentence to format.")
str_trunc1 <- function(x, n = 20, w = " ") {
    x0 <- strsplit(x, w)[[1]]
    lapply(seq(length(x0)), function(i) str_trunc0(x, i, w)) %>%
        detect(function(x) str_length(x) <= n, .dir = "backward")
}

#' @export
#' Cuts a sentence to a given number of words
#' @example str_trunc0("Hi there, I'm a sentence to format.")
str_trunc0 <- function(x, n = 5, w = " ") {
    strsplit(x, w)[[1]] %>%
        .[1:n] %>%
        paste(collapse = w)
}

str_trunc2 <- function(x, n = 20) {
    sapply(
        x,
        function(i) {
            i <- str_trim(i)
            res <- str_trunc1(i, n)
            if (is.null(res))
                res <- str_trunc1(i, n, "-")
            if (is.null(res))
                return(res)
            if (str_width(res) < str_width(i))
                paste0(res, "...")
            else
                res
        }
    ) %>% unname()
}

geneInPath <- function(x, y) {
    list.mapv(
        x[[1]] %>% as.data.frame() %>% pull(1),
        f(i) ~
            y[pull(y, 1) %in% i, ] %>%
            pull(2) %>%
            length()
    )
}

print_enrich <- function(x, path2gene = NULL, type = "enrichr", regex = NULL, pval = 0.05) {
    x <- as_tibble(x)
    if (type == "enrichr" && is.null(path2gene)) {
        stop("`path2gene` parameter must not be empty.")
    }
    if (!is.null(path2gene)) {
        x <- mutate(x, bg = geneInPath(x, path2gene))
    }
    if (type  == "enrichr") {
        x <- mutate(
            x,
            FDR = p.adjust,
            Genes = str_replace_all(geneID, "/", ";")
        )  %>% select(-geneID)
    } else if (type == "gsea") {
        x <- mutate(
            x,
            FDR = p.adjust,
            # Count = stri_split_fixed(core_enrichment, "/") %>%
            #   sapply(length),
            Count = stri_split_fixed(core_enrichment, "/") %>%
                sapply(function(x) unique(x) %>% length()),
            bg = setSize,
            Genes = stri_split_fixed(core_enrichment, "/") %>% sapply(function(x) unique(x) %>% paste0(collapse = ";"))
        )
    } else {
        x <- mutate(
            x,
            FDR = Adjusted.P.value,
            Description = Term,
            Count = str_split_fixed(Overlap, "/", 2) %>% .[, 1] %>% as.numeric(),
            bg = str_split_fixed(Overlap, "/", 2) %>% .[, 2] %>% as.numeric(),
        )
    }
    if (!is.null(regex)) {
        x <- x %>% filter(str_detect(Description, paste(regex, collapse = "|")))
    }
    x <- mutate(x, Description = to_title(Description))
    res <- x %>%
        filter(FDR <= pval) %>%
        mutate(
            # FDR = -log10(FDR) %>% round(1),
            FDR = format(FDR, digits = 3, scientific = TRUE),
            `DEG/Genes` = round(Count/bg * 100, 1)
        ) %>%
        rename(`Nb DEG` = Count, `Nb genes` = bg) %>%
        select(Description, FDR, `Nb DEG`, `Nb genes`, `DEG/Genes`, Genes, contains(c("NES", "ID")))
    if (type == "gsea") {
        res <- rename_with(
            res,
            ~ str_replace_all(., "DEG", "Enriched"),
            contains("DEG")
        ) %>%
            relocate("NES", .after = "FDR")%>%
            relocate("ID", .after = "Description")
    }
    return(res)
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
        theme_perso0(cex * 1.5) +
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
        p <- p + scale_x_continuous(labels = label_percent(1))
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
            labels = scales::label_number(accuracy = 1),
            name = title_size
        ) +
        guides(
            size = guide_legend(
                order = 1
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

plot_enrich <- function(
        x,
        n = 20,
        title = NULL,
        type = "go",
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
            str_remove_all("the ")
    }
    if (type == "gsea") {
        if (!is.null(regex)) {
            x <- filter(x, str_detect(Description, regex))
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
    } else if (type == "kegg") {
        title_size <- "# DEG"
        label_x <- "Gene ratio"
        if (!is.null(regex)) {
            x <- filter(x, str_detect(Description, regex))
        }
        df <- mutate(
            x,
            Adjusted.P.value = p.adjust,
            Term = Description %>%
                func() %>%
                to_title()
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
            x <- filter(x, str_detect(Term, regex))
        }
        df <- mutate(
            x,
            Count = str_remove_all(Overlap, "\\/.*") %>% as.numeric(),
            generatio = {
                str_split(Overlap, "/") %>%
                    sapply(function(i) as.numeric(i[1]) / as.numeric(i[2]))
            }
        )
        x_var <- "Adjusted.P.value"
    }
    df0 <- filter(df, Adjusted.P.value <= 0.05) %>%
        filter(!is.na(Term))
    if (nrow(df0) < n)
        df0 <- df
    df0 <- arrange(df0, abs(!!sym(x_var)))
    y <- "generatio"
    # if (type %in% c("gsea", "kegg")) {
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
                    str_remove_all("((ORPHA)|(WP)|(HSA)|(R-)).*") %>%
                    str_trunc2(width) %>%
                    str_trim() %>%
                    to_title()
            },
            rank = rev(row_number(!!sym(y)))
        )
    # if (type %in% c("gsea", "kegg")) {
    #   df <- mutate(df, rank = rev(row_number(!!sym(y))))
    # } else {
    df <- mutate(df, rank = row_number(!!sym(y)))
    # }
    # print(as_tibble(df) %>% select(1, 2, 4))
    if (type == "gsea") {
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


cell_paths <- c("eutrophil", "(acrophage)|(onocyte)", "endritic cell", "((natural killer)|(NK)) cell", "(T [- ]? cell)|(T-helper)|(CD[48][- ])", "B[- ]?cell", "NETosis", "Th\\d{1,2} cell")
cytok_paths0 <- c("(nterleukins?)|(IL-?\\d{1,2})", "(nterferon)|(IFN[ABG])", "(rostaglandin)|([Aa]rachidonic)|(icosa)|([Ll]eukotriene)|([Dd]ocosahexaenoic)|([Ii]cosapentaenoic)|([Ll]ipoxin)|(esolvin)")
cytok_paths <- c(cytok_paths0, "(tumor necrosis factor)|(TNF)|(NF-k)")

name_immu0 <- c(cell_paths, cytok_paths0, "STAT[ 35]", "AGE", "[Ll]upus", "(steo[cb]last)|([Bb]one)|(keletal)|(ossification)", "[Aa]rthrit", "[Gg]lucocorticoid")

name_immu <- c(name_immu0, cytok_paths, "(omplement)|([^ ]C2 )", "(oll-like)|(TLR )", "mTORC1", "(Fc gamma)|(FCG)", "etalloproteinas", "[Aa]cute") %>% unique()

name_immu2 <- c(name_immu,  "(PUMA)|(TP53)|( p53)", "([Ii]nflamm)|([Ii]mmun)", "hemokine", "mhc",  "phago((cytosis)|(some))", "leukocyte", "myeloid", "cytokine[^sis]", "granulocyte", "[Ll]ympho", "[Hh]emopo")

name_immu3 <- c(name_immu2, c("derm", "tissue", "muscl", "platelet"))

format_path <- function(x, width = 20) {
    str_remove_all(x, "\\(.*\\)") %>%
        str_remove_all("((ORPHA)|(WP)|(HSA)|(R-)).*") %>%
        # str_remove_all("(Signaling)|(By)") %>%
        # str_replace_all("Interleukin", "IL") %>%
        # str_replace(" [aA]nd ", "\\/")
        str_trunc2(width) %>%
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
    print(len0)
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
        left_join(gene_path1) %>%
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
