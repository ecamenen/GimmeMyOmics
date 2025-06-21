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
