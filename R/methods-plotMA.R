#' Plot Mean Average
#'
#' @rdname plotMA
#' @name plotMA
#'
#' @family Differential Expression Plots
#' @author Rory Kirchner, Michael Steinbaugh
#' @inheritParams AllGenerics
#'
#' @importFrom BiocGenerics plotMA
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
#' res <- examples[["res"]]
#' genes <- c("ENSMUSG00000104523", "ENSMUSG00000016918")
#'
#' # DESeqResults
#' plotMA(res, labelPoints = genes)
#'
#' # data.frame
#' df <- as.data.frame(res)
#' plotMA(df, alpha = 0.01)
NULL



# Constructors ====
#' @importFrom basejump camel
#' @importFrom dplyr filter
#' @importFrom ggplot2 aes_ annotation_logticks geom_point ggtitle guides labs
#'   scale_color_manual scale_x_log10
#' @importFrom ggrepel geom_text_repel
#' @importFrom grid arrow unit
#' @importFrom tibble rownames_to_column
.plotMA <- function(
    object,
    alpha = 0.01,
    labelPoints = NULL,
    labelColumn = "ensgene",
    pointColorScale = c("darkgray", "red"),
    labelColor = "black",
    title = TRUE) {
    results <- object %>%
        as.data.frame() %>%
        rownames_to_column("ensgene") %>%
        camel(strict = FALSE) %>%
        filter(!is.na(.data[["padj"]]))
    p <- ggplot(
        results,
        mapping = aes_(
            x = ~baseMean,
            y = ~log2FoldChange,
            color = ~padj < alpha)
    ) +
        geom_point(size = 0.8) +
        scale_x_log10() +
        annotation_logticks(sides = "b") +
        guides(color = FALSE) +
        labs(x = "mean expression across all samples",
             y = "log2 fold change")
    if (isTRUE(title)) {
        p <- p + ggtitle("ma")
    } else if (is.character(title)) {
        p <- p + ggtitle(paste("ma:", title))
    }
    if (!is.null(pointColorScale)) {
        p <- p + scale_color_manual(values = pointColorScale)
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
                force = 1,
                point.padding = unit(0.75, "lines"),
                segment.color = "gray",
                segment.size = 0.5,
                show.legend = FALSE,
                size = 4)
    }
    p
}



# Methods ====
#' @rdname plotMA
#' @importFrom ggplot2 scale_color_manual
#' @importFrom S4Vectors metadata
#' @export
setMethod(
    "plotMA",
    signature("DESeqResults"),
    function(
        object,
        labelPoints = NULL,
        labelColumn = "ensgene",
        pointColorScale = c("darkgray", "red"),
        labelColor = "black") {
        alpha <- metadata(object)[["alpha"]]
        title <- .resContrastName(object)
        .plotMA(
            object,
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
setMethod(
    "plotMA",
    signature("data.frame"),
    .plotMA)
