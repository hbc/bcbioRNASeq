#' Plot row standard deviations versus row means
#'
#' Modified wrapper function of [vsn::meanSdPlot()] that plots [log2()],
#' [rlog()], and [varianceStabilizingTransformation()] count data.
#'
#'
#' @rdname qc_plots
#' @param dds [DESeqDataSet].
#' @param rld [DESeqTransform] resulting from [rlog()] applied to
#'   [DESeqDataSet].
#' @param vsd [DESeqTransform] resulting from
#'   [varianceStabilizingTransformation()] applied to [DESeqDataSet].
#'
#' @export
plot_mean_sd <- function(dds, rld = NULL, vsd = NULL) {
    nonzero <- dds %>% counts %>% rowSums %>% `>`(0)

    # log2
    dds %>%
        counts(normalized = TRUE) %>%
        .[nonzero, ] %>%
        `+`(1) %>%
        log2 %>%
        meanSdPlot

    # rlog transformation
    # `rld` = RLog Data
    if (is.null(rld)) {
        rld <- rlog(dds)
    }
    rld %>%
        assay %>%
        .[nonzero, ] %>%
        meanSdPlot

    # Variance stabilizing transformation (VST)
    # `vsd` = VSt Data
    if (is.null(vsd)) {
        vsd <- varianceStabilizingTransformation(dds)
    }
    vsd %>%
        assay %>%
        .[nonzero, ] %>%
        meanSdPlot
}
