#' Plot Row Standard Deviations vs. Row Means
#'
#' [vsn::meanSdPlot()] wrapper that plots [log2()], [rlog()], and
#' [varianceStabilizingTransformation()] normalized counts.
#'
#' @rdname plotMeanSD
#' @name plotMeanSD
#' @family Differential Expression Utilities
#' @author Michael Steinbaugh, Lorena Patano
#'
#' @inheritParams AllGenerics
#'
#' @param orientation Orientation to use for plot grid, either `horizontal` or
#'   `vertical` (default).
#' @param showLegend Include the color bar legend. This is typically not that
#'   informative and is disabled by default, to improve the plot appearance.
#'
#' @return [ggplot] grid.
#'
#' @examples
#' data(bcb, dds)
#'
#' # bcbioRNASeq
#' plotMeanSD(bcb)
#'
#' # DESeqDataSet
#' plotMeanSD(dds)
NULL



# Constructors ====
.plotMeanSD <- function(
    raw,
    normalized,
    rlog,
    vst,
    orientation = "vertical",
    showLegend = FALSE) {
    xlab <- "rank (mean)"
    nonzero <- raw %>%
        rowSums() %>%
        `>`(0)
    gglog2 <- normalized %>%
        .[nonzero, , drop = FALSE] %>%
        `+`(1) %>%
        log2() %>%
        meanSdPlot(plot = FALSE) %>%
        .[["gg"]] +
        ggtitle("log2") +
        xlab(xlab)
    ggrlog <- rlog %>%
        .[nonzero, , drop = FALSE] %>%
        meanSdPlot(plot = FALSE) %>%
        .[["gg"]] +
        ggtitle("rlog") +
        xlab(xlab) +
        theme(legend.position = "none")
    ggvst <- vst %>%
        .[nonzero, , drop = FALSE] %>%
        meanSdPlot(plot = FALSE) %>%
        .[["gg"]] +
        ggtitle("vst") +
        xlab(xlab) +
        theme(legend.position = "none")

    # Remove the plot (color) legend, if desired
    if (!isTRUE(showLegend)) {
        gglog2 <- gglog2 +
            theme(legend.position = "none")
        ggrlog <- ggrlog +
            theme(legend.position = "none")
        ggvst <- ggvst +
            theme(legend.position = "none")
    }

    # Return either horizontal or vertical
    if (orientation == "horizontal") {
        ncol <- 3
        nrow <- 1
    } else if (orientation == "vertical") {
        ncol <- 1
        nrow <- 3
    }

    plot_grid(
        gglog2,
        ggrlog,
        ggvst,
        labels = "AUTO",
        ncol = ncol,
        nrow = nrow)
}


# Methods ====
#' @rdname plotMeanSD
#' @export
setMethod(
    "plotMeanSD",
    signature("bcbioRNASeqANY"),
    function(
        object,
        orientation = "vertical",
        showLegend = FALSE) {
        .plotMeanSD(
            raw = counts(object, normalized = FALSE),
            normalized = counts(object, normalized = TRUE),
            rlog = counts(object, normalized = "rlog"),
            vst = counts(object, normalized = "vst"))
    })



#' @rdname plotMeanSD
#' @export
setMethod(
    "plotMeanSD",
    signature("DESeqDataSet"),
    function(
        object,
        orientation = "vertical",
        showLegend = FALSE) {
        .plotMeanSD(
            raw = counts(object, normalized = FALSE),
            normalized = counts(object, normalized = TRUE),
            rlog = assay(rlog(object)),
            vst = assay(varianceStabilizingTransformation(object)))
    })
