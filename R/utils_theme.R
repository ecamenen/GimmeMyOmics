theme_custom <- function(
        cex = 1,
        cex_main = 17 * cex,
        cex_sub = 15 * cex,
        cex_axis = 15 * cex) {
    theme_minimal() +
        theme(
            axis.text = element_text(size = 13 * cex, color = "gray50"),
            axis.title = element_text(face = "bold.italic", size = cex_axis),
            strip.text = element_text(
                size = cex_main,
                face = "bold",
                hjust = 0.5,
                margin = margin(0.5, 0.5, 0.5, 0.5)
            ),
            plot.title = element_text(face = "bold", size = cex_main, hjust = 0.5),
            legend.title = element_text(face = "italic", size = cex_sub),
            legend.text = element_text(colour = "black", size = 10 * cex)
        )
}

kable0 <- function(x, align = "c", color = "#a9a9a9") {
    x %>%
        kbl(escape = FALSE, align = c("l", rep(align, ncol(x) -1))) %>%
        kable_minimal(full_width = FALSE) %>%
        column_spec(1, bold = TRUE, color = color)
}

kable2 <- function(x, align = "c", color = "#a9a9a9")
    x %>%
    set_colnames(colnames(.) %>% str_clean()) %>%
    kable0(align, color)

list_cbind0 <- function(x) {
    row_names <- lapply(x, rownames)
    x <- list.map(x, f(i) ~ as.data.frame(i))
    tmp <- sapply(x, nrow)
    for (i in seq(length(tmp))) {
        add <- NULL
        add0 <- NULL
        if(i > 1) {
            add <- matrix(nrow = sum(tmp[1:(i-1)]), ncol = ncol(x[[i]]))
            if (is.data.frame(x[[i]])) {
                add <- as.data.frame(add)
            }
            colnames(add) <- colnames(x[[i]])
        }
        if(i < length(tmp)) {
            add0 <- matrix(nrow = sum(tmp[(i+1):length(tmp)]), ncol = ncol(x[[i]]))
            if (is.data.frame(x[[i]])) {
                add0 <- as.data.frame(add0)
            }
            colnames(add0) <- colnames(x[[i]])
        }
        x[[i]] <- rbind(add, x[[i]], add0)
        rownames(x[[i]]) <- unlist(row_names)
    }
    list.cbind(x)
}
