#' Plot Row Standard Deviations vs. Row Means
#'
#' [vsn::meanSdPlot()] wrapper that plots [log2()], [rlog()], and
#' [varianceStabilizingTransformation()] normalized counts.
#'
#' @rdname plotMeanSD
#' @name plotMeanSD
#'
#' @return [ggplot] grid.
#' @export
#'
#' @examples
#' data(bcb, dds)
#'
#' # bcbioRNADataSet
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
    nonzero <- raw %>%
        rowSums %>%
        `>`(0L)
    gglog2 <- normalized %>%
        .[nonzero, ] %>%
        `+`(1L) %>%
        log2 %>%
        meanSdPlot(plot = FALSE)
    ggrlog <- rlog %>%
        .[nonzero, ] %>%
        meanSdPlot(plot = FALSE)
    ggvst <- vst %>%
        .[nonzero, ] %>%
        meanSdPlot(plot = FALSE)
    plot_grid(
        gglog2[["gg"]]  +
            ggtitle("log2") +
            theme(legend.position = "none"),
        ggrlog[["gg"]] +
            ggtitle("rlog") +
            theme(legend.position = "none"),
        ggvst[["gg"]] +
            ggtitle("variance stabilizing transformation") +
            theme(legend.position = "none"),
        labels = "auto",
        nrow = 3L)
}


# Methods ====
#' @rdname plotMeanSD
#' @export
setMethod("plotMeanSD", "bcbioRNADataSet", function(object) {
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
        rlog = rlog(object) %>% assay,
        vst = varianceStabilizingTransformation(object) %>% assay)
})
