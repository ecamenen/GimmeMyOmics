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

kable0 <- function(x, align = "c", color = "#a9a9a9", func = str_clean) {
    x %>%
        set_colnames(colnames(.) %>% func()) %>%
        kbl(escape = FALSE, align = c("l", rep(align, ncol(x) -1))) %>%
        kable_minimal(full_width = FALSE) %>%
        column_spec(1, bold = TRUE, color = color)
}
