#' @name plotMeanSD
#' @inherit bioverbs::plotMeanSD
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Lorena Patano
#'
#' @details
#' [vsn::meanSdPlot()] wrapper that plots count transformations on a log2 scale.
#'
#' - **DESeq2 log2**: log2 normalized counts.
#' - **DESeq2 rlog**: regularized log transformation.
#' - **DESeq2 vst**: variance stabilizing transformation.
#' - **edgeR log2 tmm**: log2 trimmed mean of M-values transformation.
#'
#' @seealso
#' - [vsn::meanSdPlot()].
#' - [DESeq2::DESeq()].
#' - [DESeq2::rlog()].
#' - [DESeq2::varianceStabilizingTransformation()].
#' - [tmm()].
#'
#' @inheritParams general
#' @param ... Additional arguments.
#'
#' @return `ggplot` grid.
#'
#' @examples
#' # bcbioRNASeq ====
#' plotMeanSD(bcb_small)
#'
#' # DESeqDataSet ====
#' plotMeanSD(dds_small)
NULL



#' @rdname plotMeanSD
#' @name plotMeanSD
#' @importFrom bioverbs plotMeanSD
#' @usage plotMeanSD(object, ...)
#' @export
NULL



.plotMeanSD <- function(
    raw,
    normalized,
    rlog,
    vst,
    legend = FALSE
) {
    assert_is_matrix(raw)
    assert_is_matrix(normalized)
    # rlog and vst are optional (transformationLimit)
    assert_is_a_bool(legend)

    xlab <- "rank (mean)"
    nonzero <- rowSums(raw) > 0L

    # DESeq2 log2 normalized
    gglog2 <- normalized %>%
        .[nonzero, , drop = FALSE] %>%
        `+`(1L) %>%
        log2() %>%
        meanSdPlot(plot = FALSE) %>%
        .[["gg"]] +
        ggtitle("DESeq2 log2") +
        xlab(xlab)

    # DESeq2 regularized log
    if (is.matrix(rlog)) {
        ggrlog <- rlog %>%
            .[nonzero, , drop = FALSE] %>%
            meanSdPlot(plot = FALSE) %>%
            .[["gg"]] +
            ggtitle("DESeq2 rlog") +
            xlab(xlab)
    } else {
        message("Skipping regularized log")
        ggrlog <- NULL
    }

    # DESeq2 variance stabilizing transformation
    if (is.matrix(vst)) {
        ggvst <- vst %>%
            .[nonzero, , drop = FALSE] %>%
            meanSdPlot(plot = FALSE) %>%
            .[["gg"]] +
            ggtitle("DESeq2 vst") +
            xlab(xlab)
    } else {
        message("Skipping variance stabilizing transformation")
        ggvst <- NULL
    }

    # edgeR log2 tmm
    tmm <- suppressMessages(tmm(raw))
    ggtmm <- tmm %>%
        .[nonzero, , drop = FALSE] %>%
        `+`(1L) %>%
        log2() %>%
        meanSdPlot(plot = FALSE) %>%
        .[["gg"]] +
        ggtitle("edgeR log2 tmm") +
        xlab(xlab)

    plotlist <- list(
        log2 = gglog2,
        rlog = ggrlog,
        vst = ggvst,
        tmm = ggtmm
    )
    # Remove NULL assays if present (e.g. rlog, vst)
    plotlist <- Filter(Negate(is.null), plotlist)

    # Remove the plot (color) legend, if desired
    if (!isTRUE(legend)) {
        plotlist <- lapply(plotlist, function(p) {
            p <- p + theme(legend.position = "none")
        })
    }

    plot_grid(plotlist = plotlist)
}



# Require that the DESeq2 transformations are slotted.
# If `transformationLimit` was applied, this function will error.

plotMeanSD.bcbioRNASeq <-  # nolint
    function(
        object,
        legend = getOption("bcbio.legend", FALSE)
    ) {
        .plotMeanSD(
            raw = counts(object, normalized = FALSE),
            normalized = counts(object, normalized = TRUE),
            rlog = assays(object)[["rlog"]],
            vst = assays(object)[["vst"]],
            legend = legend
        )
    }



#' @rdname plotMeanSD
#' @export
setMethod(
    f = "plotMeanSD",
    signature = signature("bcbioRNASeq"),
    definition = plotMeanSD.bcbioRNASeq
)



plotMeanSD.DESeqDataSet <-  # nolint
    function(
        object,
        legend = getOption("bcbio.legend", FALSE)
    ) {
        .plotMeanSD(
            raw = counts(object, normalized = FALSE),
            normalized = counts(object, normalized = TRUE),
            rlog = assay(rlog(object)),
            vst = assay(varianceStabilizingTransformation(object)),
            legend = legend
        )
    }



#' @rdname plotMeanSD
#' @export
setMethod(
    f = "plotMeanSD",
    signature = signature("DESeqDataSet"),
    definition = plotMeanSD.DESeqDataSet
)
