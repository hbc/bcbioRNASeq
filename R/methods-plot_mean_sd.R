#' Plot row standard deviations versus row means
#'
#' [vsn::meanSdPlot()] wrapper that plots [log2()], [rlog()], and
#' [varianceStabilizingTransformation()] normalized counts.
#'
#' @rdname plot_mean_sd
#' @author Michael Steinbaugh, Lorena Patano
#' @family DESeq2 Utilities
#'
#' @return [ggplot] grid.
#' @export
#'
#' @examples
#' data(bcb)
#' plot_mean_sd(bcb)
setMethod("plot_mean_sd", "bcbioRNADataSet", function(object) {
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
    plot_grid(gglog2[["gg"]]  + theme(legend.position = "none"),
              ggrlog[["gg"]] + theme(legend.position = "none"),
              ggvsd[["gg"]] + theme(legend.position = "none"),
              nrow = 3L)
})
