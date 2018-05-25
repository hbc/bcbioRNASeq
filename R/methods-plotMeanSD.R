#' Plot Row Standard Deviations vs. Row Means
#'
#' [vsn::meanSdPlot()] wrapper that plots [log2()], [rlog()], and
#' [varianceStabilizingTransformation()] normalized counts.
#'
#' @name plotMeanSD
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Lorena Patano
#'
#' @inheritParams general
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



# Constructors =================================================================
.plotMeanSD <- function(
    raw,
    normalized,
    rlog,
    vst,
    legend = FALSE
) {
    assert_is_matrix(raw)
    assert_is_matrix(normalized)
    assert_is_matrix(rlog)
    assert_is_matrix(vst)
    assert_is_a_bool(legend)

    xlab <- "rank (mean)"
    nonzero <- rowSums(raw) > 0L

    gglog2 <- normalized %>%
        .[nonzero, , drop = FALSE] %>%
        `+`(1L) %>%
        log2() %>%
        meanSdPlot(plot = FALSE) %>%
        .[["gg"]] +
        ggtitle("log2") +
        xlab(xlab)
    tmm <- suppressMessages(tmm(raw))
    ggtmm <- tmm %>%
        .[nonzero, , drop = FALSE] %>%
        `+`(1L) %>%
        log2() %>%
        meanSdPlot(plot = FALSE) %>%
        .[["gg"]] +
        ggtitle("log2 tmm") +
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
        ggtitle("vst") +
        xlab(xlab)

    plotlist <- list(
        log2 = gglog2,
        tmm = ggtmm,
        rlog = ggrlog,
        vst = ggvst
    )

    # Remove the plot (color) legend, if desired
    if (!isTRUE(legend)) {
        plotlist <- lapply(plotlist, function(p) {
            p <- p + theme(legend.position = "none")
        })
    }

    plot_grid(plotlist = plotlist, labels = "AUTO")
}



# Methods ======================================================================
# Require that the DESeq2 transformations are slotted.
# If `transformationLimit` was applied, this function will error.
#' @rdname plotMeanSD
#' @export
setMethod(
    "plotMeanSD",
    signature("bcbioRNASeq"),
    function(
        object,
        legend = FALSE
    ) {
        .plotMeanSD(
            raw = counts(object, normalized = FALSE),
            normalized = counts(object, normalized = TRUE),
            rlog = assays(object)[["rlog"]],
            vst = assays(object)[["vst"]],
            legend = legend
        )
    }
)



#' @rdname plotMeanSD
#' @export
setMethod(
    "plotMeanSD",
    signature("DESeqDataSet"),
    function(
        object,
        legend = FALSE
    ) {
        .plotMeanSD(
            raw = counts(object, normalized = FALSE),
            normalized = counts(object, normalized = TRUE),
            rlog = assay(rlog(object)),
            vst = assay(varianceStabilizingTransformation(object)),
            legend = legend
        )
    }
)
