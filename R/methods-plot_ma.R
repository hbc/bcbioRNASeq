#' MA-plot from base means and log fold changes
#'
#' @rdname plot_ma
#' @author Michael Steinbaugh, Rory Kirchner
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
.plot_ma <- function(res_df, title, label_points, label_column) {
    res_df <- filter(res_df, !is.na(.data[["padj"]]))
    p <- ggplot(res_df,
                aes_(x = ~baseMean,
                     y = ~log2FoldChange,
                     color = ~padj < 0.05)) +
        geom_point(size = 0.8) +
        scale_x_log10(
            breaks = trans_breaks("log10", function(x) 10L ^ x),
            labels = trans_format("log10", math_format(10L ^ .x))) +  # nolint
        annotation_logticks(sides = "b") +
        xlab("mean expression across all samples") +
        ylab(expression(log[2]*" fold change")) +  # nolint
        scale_color_manual(values = c("black", "red", "green")) +
        guides(color = FALSE)
    if (!is.null(title)) {
        p <- p + ggtitle(title)
    }
    if (!is.null(label_points)) {
        labels <- res_df[res_df[[label_column]] %in% label_points, ]
        p <- p + geom_text(data = labels,
                           aes_string("baseMean", "log2FoldChange",
                                      label = label_column), size = 3L)
    }
    p
}



#' @rdname plot_ma
#' @export
setMethod(
    "plot_ma",
    "DESeqResults",
    function(object,
             title = NULL,
             label_points = NULL,
             label_column = "symbol") {
        res %>%
            as.data.frame %>%
            .plot_ma(title = title,
                     label_points = label_points,
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
