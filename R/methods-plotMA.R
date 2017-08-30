#' Plot Mean Average
#'
#' @rdname plotMA
#' @name plotMA
#'
#' @family Differential Expression Plots
#' @author Rory Kirchner, Michael Steinbaugh
#' @inheritParams AllGenerics
#'
#' @param alpha Alpha level cutoff (Adjusted P value).
#' @param labelPoints *Optional*. Label these particular points.
#' @param labelColumn Match `labelPoints` argument to this column in the
#'   results.
#' @param pointColorScale Point color scale. See
#'   [ggplot2::scale_color_manual()] for more information.
#' @param labelColor Text label color.
#' @param title Plot title.
#'
#' @return [ggplot].
#'
#' @examples
#' data(res)
#' genes <- c("ENSMUSG00000104523", "ENSMUSG00000016918")
#'
#' # DESeqResults
#' plotMA(res, labelPoints = genes)
#'
#' # data.frame
#' df <- as.data.frame(res)
#' plotMA(df, alpha = 0.1, labelPoints = genes)
NULL



# Constructors ====
.plotMA <- function(
    object,
    alpha = 0.01,
    labelPoints = NULL,
    labelColumn = "rowname",
    pointColorScale = c("darkgrey", "red", "green"),
    labelColor = "black",
    title = TRUE) {
    results <- object %>%
        as("tibble") %>%
        camel %>%
        .[!is.na(.[["padj"]]), ]
    p <- ggplot(results,
                aes_(x = ~baseMean,
                     y = ~log2FoldChange,
                     color = ~padj < alpha)) +
        geom_point(size = 0.8) +
        scale_x_log10() +
        annotation_logticks(sides = "b") +
        scale_color_manual(values = pointColorScale) +
        guides(color = FALSE) +
        labs(x = "mean expression across all samples",
             y = "log2 fold change")
    if (isTRUE(title)) {
        p <- p +
            ggtitle("ma")
    } else if (is.character(title)) {
        p <- p +
            ggtitle(paste("ma:", title))
    }
    if (!is.null(labelPoints)) {
        labels <- results %>%
            .[.[[labelColumn]] %in% labelPoints, , drop = FALSE]
        p <- p +
            geom_text_repel(
                data = labels,
                aes_(x = ~baseMean,
                     y = ~log2FoldChange,
                     label = as.name(labelColumn)),
                arrow = arrow(length = unit(0.01, "npc")),
                box.padding = unit(0.5, "lines"),
                color = labelColor,
                fontface = "bold",
                force = 1L,
                point.padding = unit(0.75, "lines"),
                segment.color = "gray",
                segment.size = 0.5,
                show.legend = FALSE,
                size = 4L)
    }
    p
}



# Methods ====
#' @rdname plotMA
#' @export
setMethod("plotMA", "DESeqResults", function(
    object,
    labelPoints = NULL,
    labelColumn = "rowname",
    pointColorScale = c("darkgrey", "red", "green"),
    labelColor = "black") {
    alpha <- metadata(object)[["alpha"]]
    title <- .resContrastName(object)
    .plotMA(object,
            labelPoints = labelPoints,
            labelColumn = labelColumn,
            pointColorScale = pointColorScale,
            labelColor = labelColor,
            # Automatic
            alpha = alpha,
            title = title)
})



#' @rdname plotMA
#' @export
setMethod("plotMA", "data.frame", .plotMA)
