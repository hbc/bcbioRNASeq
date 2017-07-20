#' Plot Mean Average
#'
#' @rdname plot_ma
#' @author Michael Steinbaugh, Rory Kirchner
#' @family Differential Expression Plots
#'
#' @param object Object.
#' @param title *Optional*. Plot title.
#' @param label_points *Optional*. Label these particular points.
#' @param label_column Match `label_points` argument to this column in the
#'   results.
#' @param point_color_scale Point color scale. See [ggplot2::color_manual] for
#'   more information.
#' @param label_color Text label color.
#'
#' @return [ggplot].
#'
#' @examples
#' data(bcb)
#' dds <- DESeqDataSetFromTximport(
#'     txi = txi(bcb),
#'     colData = colData(bcb),
#'     design = formula(~group)) %>%
#'     DESeq
#' res <- results(dds)
#'
#' genes <- c("ENSMUSG00000104523", "ENSMUSG00000016918")
#'
#' # DESeqResults
#' plot_ma(res, label_points = genes)
#'
#' # data.frame
#' res_df <- as.data.frame(res)
#' plot_ma(res_df, label_points = genes)



#' @rdname plot_ma
#' @usage NULL
.plot_ma <- function(
    res_df,
    alpha = 0.05,
    label_points = NULL,
    label_column = "rowname",
    point_color_scale = c("darkgrey", "red", "green"),
    label_color = "black",
    title = TRUE) {
    res_tbl <- res_df %>%
        as("tibble") %>%
        snake %>%
        filter(!is.na(.data[["padj"]]))
    p <- ggplot(res_tbl,
                aes_(x = ~base_mean,
                     y = ~log2_fold_change,
                     color = ~padj < alpha)) +
        geom_point(size = 0.8) +
        scale_x_log10() +
        annotation_logticks(sides = "b") +
        scale_color_manual(values = point_color_scale) +
        guides(color = FALSE) +
        labs(x = "mean expression across all samples",
             y = "log2 fold change")
    if (isTRUE(title)) {
        p <- p + ggtitle("mean average")
    } else if (is.character(title)) {
        p <- p + ggtitle(title)
    }
    if (!is.null(label_points)) {
        labels <- res_tbl %>%
            filter(.data[[label_column]] %in% !!syms(label_points))
        p <- p +
            geom_text_repel(
                data = labels,
                aes_(x = ~base_mean,
                     y = ~log2_fold_change,
                     label = as.name(label_column)),
                color = label_color,
                size = 4L)
    }
    p
}



#' @rdname plot_ma
#' @export
setMethod("plot_ma", "DESeqResults", function(object, ...) {
    res %>% as.data.frame %>% .plot_ma(...)
})



#' @rdname plot_ma
#' @export
setMethod("plot_ma", "data.frame", function(object, ...) {
    .plot_ma(object, ...)
})
