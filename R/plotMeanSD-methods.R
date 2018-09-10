#' Plot Row Standard Deviations vs. Row Means
#'
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
#' @name plotMeanSD
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Lorena Patano
#' @export
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
#' dds_small <- as(deseq_small, "DESeqDataSet")
#' plotMeanSD(dds_small)
NULL



.plotMeanSD <- function(
    raw,
    normalized,
    vst = NULL,
    rlog = NULL,
    legend = FALSE
) {
    assert_is_matrix(raw)
    assert_is_matrix(normalized)
    assert_are_identical(dimnames(normalized), dimnames(raw))
    assert_is_any_of(vst, c("matrix", "NULL"))
    if (is.matrix(vst)) {
        assert_are_identical(dimnames(vst), dimnames(raw))
    }
    assert_is_any_of(rlog, c("matrix", "NULL"))
    if (is.matrix(rlog)) {
        assert_are_identical(dimnames(rlog), dimnames(raw))
    }
    assert_is_a_bool(legend)

    xlab <- "rank (mean)"
    nonzero <- rowSums(raw) > 0L

    # DESeq2 log2 normalized.
    gglog2 <- normalized %>%
        .[nonzero, , drop = FALSE] %>%
        `+`(1L) %>%
        log2() %>%
        meanSdPlot(plot = FALSE) %>%
        .[["gg"]] +
        ggtitle("DESeq2 log2") +
        xlab(xlab)

    # DESeq2 regularized log.
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

    # DESeq2 variance stabilizing transformation.
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

    # edgeR log2 tmm.
    ggtmm <- raw %>%
        tmm() %>%
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
    # Remove NULL assays if present (e.g. rlog, vst).
    plotlist <- Filter(Negate(is.null), plotlist)

    # Remove the plot (color) legend, if desired.
    if (!isTRUE(legend)) {
        plotlist <- lapply(plotlist, function(p) {
            p <- p + theme(legend.position = "none")
        })
    }

    plot_grid(plotlist = plotlist)
}



# Require that the DESeq2 transformations are slotted.
# If `transformationLimit` was applied, this function will error.
#' @rdname plotMeanSD
#' @export
setMethod(
    "plotMeanSD",
    signature("bcbioRNASeq"),
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
)



#' @rdname plotMeanSD
#' @export
setMethod(
    "plotMeanSD",
    signature("DESeqDataSet"),
    function(
        object,
        vst = TRUE,
        rlog = FALSE,
        legend = getOption("bcbio.legend", FALSE)
    ) {
        validObject(object)
        assert_is_a_bool(vst)
        if (isTRUE(vst)) {
            vst <- assay(varianceStabilizingTransformation(object))
        } else {
            vst <- NULL
        }
        assert_is_a_bool(rlog)
        if (isTRUE(rlog)) {
            rlog <- assay(rlog(object))
        } else {
            rlog <- NULL
        }
        .plotMeanSD(
            raw = counts(object, normalized = FALSE),
            normalized = counts(object, normalized = TRUE),
            vst = vst,
            rlog = rlog,
            legend = legend
        )
    }
)
