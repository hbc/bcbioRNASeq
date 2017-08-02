#' Plot Row Standard Deviations vs. Row Means
#'
#' [vsn::meanSdPlot()] wrapper that plots [log2()], [rlog()], and
#' [varianceStabilizingTransformation()] normalized counts.
#'
#' @rdname plotMeanSD
#' @author Michael Steinbaugh, Lorena Patano
#' @family DESeq2 Utilities
#'
#' @return [ggplot] grid.
#' @export
#'
#' @examples
#' data(bcb)
#' plotMeanSD(bcb)
setMethod("plotMeanSD", "bcbioRNADataSet", function(object) {
    nonzero <- counts(object, normalized = FALSE) %>%
        rowSums %>%
        `>`(0L)
    gglog2 <- counts(object, normalized = TRUE) %>%
        .[nonzero, ] %>%
        `+`(1L) %>%
        log2 %>%
        meanSdPlot(plot = FALSE)
    ggrlog <- counts(object, normalized = "rlog") %>%
        .[nonzero, ] %>%
        meanSdPlot(plot = FALSE)
    ggvsd <- counts(object, normalized = "vst") %>%
        .[nonzero, ] %>%
        meanSdPlot(plot = FALSE)
    plot_grid(
        gglog2[["gg"]]  +
            ggtitle("log2") +
            theme(legend.position = "none"),
        ggrlog[["gg"]] +
            ggtitle("rlog") +
            theme(legend.position = "none"),
        ggvsd[["gg"]] +
            ggtitle("variance stabilizing transformation") +
            theme(legend.position = "none"),
        labels = "auto",
        nrow = 3L)
})
