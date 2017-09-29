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
    vst) {
    xlab <- "rank (mean)"
    nonzero <- raw %>%
        rowSums %>%
        `>`(0)
    gglog2 <- normalized %>%
        .[nonzero, , drop = FALSE] %>%
        `+`(1) %>%
        log2 %>%
        meanSdPlot(plot = FALSE) %>%
        .[["gg"]] +
        ggtitle("log2") +
        xlab(xlab)
    ggrlog <- rlog %>%
        .[nonzero, , drop = FALSE] %>%
        meanSdPlot(plot = FALSE) %>%
        .[["gg"]] +
        ggtitle("rlog") +
        xlab(xlab)
    ggvst <- vst %>%
        .[nonzero, , drop = FALSE] %>%
        meanSdPlot(plot = FALSE) %>%
        .[["gg"]] +
        ggtitle("variance stabilizing transformation") +
        xlab(xlab)
    plot_grid(
        gglog2,
        ggrlog,
        ggvst,
        labels = "auto",
        nrow = 3)
}


# Methods ====
#' @rdname plotMeanSD
#' @export
setMethod("plotMeanSD", "bcbioRNASeqANY", function(object) {
    .plotMeanSD(
        raw = counts(object, normalized = FALSE),
        normalized = counts(object, normalized = TRUE),
        rlog = counts(object, normalized = "rlog"),
        vst = counts(object, normalized = "vst"))
})



#' @rdname plotMeanSD
#' @export
setMethod("plotMeanSD", "DESeqDataSet", function(object) {
    .plotMeanSD(
        raw = counts(object, normalized = FALSE),
        normalized = counts(object, normalized = TRUE),
        rlog = assay(rlog(object)),
        vst = assay(varianceStabilizingTransformation(object)))
})
