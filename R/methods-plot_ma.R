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
#' plot_ma(res)



#' @rdname plot_ma
#' @usage NULL
.plot_ma <- function(
    res_df,
    alpha = 0.05,
    label_points = NULL,
    label_column = "rowname") {
    res_tbl <- res_df %>%
        as("tibble") %>%
        snake %>%
        filter(!is.na(.data[["padj"]]))
    p <- ggplot(res_tbl,
                aes_(x = ~base_mean,
                     y = ~log2_fold_change,
                     color = ~padj < !!alpha)) +
        geom_point(size = 0.8) +
        scale_x_log10() +
        annotation_logticks(sides = "b") +
        scale_color_manual(values = c("black", "red", "green")) +
        guides(color = FALSE) +
        labs(title = "mean average",
             x = "mean expression across all samples",
             y = "log2 fold change")
    if (!is.null(label_points)) {
        labels <- res_tbl %>%
            .[.[[label_column]] %in% label_points, ]
        p <- p +
            geom_text_repel(
                data = labels,
                aes_string("base_mean",
                           "log2_fold_change",
                           label = label_column),
                size = 3L)
    }
    p
}



#' @rdname plot_ma
#' @export
setMethod(
    "plot_ma",
    "DESeqResults",
    function(object, label_points = NULL) {
        res %>%
            as.data.frame %>%
            # FIXME Add left_join for `label_column` = symbol here
            # Automatically label the top 30 if TRUE?
            .plot_ma(label_points = label_points,
                     label_column = label_column)
    })



#' @rdname plot_ma
#' @export
setMethod(
    "plot_ma",
    "data.frame",
    function(object,
             title = NULL,
             label_points = NULL,
             label_column = "symbol") {
        .plot_ma
    })
