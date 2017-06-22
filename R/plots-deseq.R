#' Plot row standard deviations versus row means
#'
#' Modified wrapper function of [vsn::meanSdPlot()] that plots [log2()],
#' [rlog()], and [varianceStabilizingTransformation()] count data.
#'
#' @param dds [DESeqDataSet].
#' @param rld [DESeqTransform] resulting from [rlog()] applied to
#'   [DESeqDataSet].
#' @param vsd [DESeqTransform] resulting from
#'   [varianceStabilizingTransformation()] applied to [DESeqDataSet].
#'
#' @export
plot_mean_sd <- function(dds, rld = NULL, vsd = NULL) {
    nonzero <- dds %>%
        counts %>%
        rowSums %>%
        `>`(0L)

    # log2 ----
    gglog2 <- dds %>%
        counts(normalized = TRUE) %>%
        .[nonzero, ] %>%
        `+`(1L) %>%
        log2 %>%
        meanSdPlot(plot = FALSE)

    # rlog transformation ----
    # `rld` = RLog Data
    if (is.null(rld)) {
        rld <- rlog(dds)
    }
    ggrlog <- rld %>%
        assay %>%
        .[nonzero, ] %>%
        meanSdPlot(plot = FALSE)

    # Variance stabilizing transformation (VST) ----
    # `vsd` = VSt Data
    if (is.null(vsd)) {
        vsd <- varianceStabilizingTransformation(dds)
    }
    ggvsd <- vsd %>%
        assay %>%
        .[nonzero, ] %>%
        meanSdPlot(plot = FALSE)
    plot_grid(gglog2[["gg"]]  + theme(legend.position = "none"),
              ggrlog[["gg"]] + theme(legend.position = "none"),
              ggvsd[["gg"]] + theme(legend.position = "none"), nrow = 1L)
}
