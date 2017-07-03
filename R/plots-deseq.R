#' Plot dispersion estimates
#'
#' [DESeq2::plotDispEsts()] wrapper that supports a [bcbioRNADataSet].
#'
#' @author Michael Steinbaugh
#'
#' @param bcb [bcbioRNADataSet].
#'
#' @return [ggplot].
#' @export
plot_dispersion <- function(bcb) {
    bcbio("DESeqDataSet") %>%
        plotDispEsts
}



#' Plot row standard deviations versus row means
#'
#' [vsn::meanSdPlot()] wrapper that plots [log2()], [rlog()], and
#' [varianceStabilizingTransformation()] normalized counts.
#'
#' @author Michael Steinbaugh, Lorena Patano
#'
#' @param bcb [bcbioRNADataSet].
#'
#' @return ggplot grid.
#' @export
plot_mean_sd <- function(bcb) {
    nonzero <- counts(bcb, normalized = FALSE) %>%
        rowSums %>%
        `>`(0L)
    gglog2 <- counts(bcb, normalized = TRUE) %>%
        .[nonzero, ] %>%
        `+`(1L) %>%
        log2 %>%
        meanSdPlot(plot = FALSE)
    ggrlog <- counts(bcb, normalized = "rlog") %>%
        .[nonzero, ] %>%
        meanSdPlot(plot = FALSE)
    ggvsd <- counts(bcb, normalized = "vst") %>%
        .[nonzero, ] %>%
        meanSdPlot(plot = FALSE)
    plot_grid(gglog2[["gg"]]  + theme(legend.position = "none"),
              ggrlog[["gg"]] + theme(legend.position = "none"),
              ggvsd[["gg"]] + theme(legend.position = "none"), nrow = 1L)
}
