#' @name plotMeanSD
#' @author Michael Steinbaugh, Lorena Patano
#' @inherit bioverbs::plotMeanSD
#' @inheritParams basejump::params
#' @inheritParams params
#' @inheritParams bcbioRNASeq
#'
#' @details
#' `vsn::meanSdPlot` wrapper that plots count transformations on a log2 scale.
#'
#' - DESeq2 log2: log2 size factor-adjusted normalized counts.
#' - DESeq2 rlog: **r**egularized **log** transformation.
#' - DESeq2 VST: **v**ariance-**s**tabilizing **t**ransformation.
#' - edgeR log2 TMM: log2 **t**rimmed **m**ean of **M**-values transformation.
#'
#' @seealso
#' - `vsn::meanSdPlot()`.
#' - `DESeq2::DESeq()`.
#' - `DESeq2::rlog()`.
#' - `DESeq2::varianceStabilizingTransformation()`.
#' - `edgeR::calcNormFactors()`.
#'
#' @examples
#' data(bcb)
#' plotMeanSD(bcb)
NULL



#' @importFrom bioverbs plotMeanSD
#' @aliases NULL
#' @export
bioverbs::plotMeanSD



# Match the vst, rlog conventions to `bcbioRNASeq` generator.
.plotMeanSD <- function(
    raw,
    normalized,
    vst = NULL,
    rlog = NULL,
    legend = FALSE
) {
    assert(
        is.matrix(raw),
        is.matrix(normalized),
        identical(dimnames(normalized), dimnames(raw)),
        isAny(vst, classes = c("matrix", "NULL")),
        isAny(rlog, classes = c("matrix", "NULL")),
        isFlag(legend)
    )
    if (is.matrix(vst)) {
        assert(identical(dimnames(vst), dimnames(raw)))
    }
    if (is.matrix(rlog)) {
        assert(identical(dimnames(rlog), dimnames(raw)))
    }

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
        message(paste(
            "Regularized log (rlog) was not calculated.",
            "Skipping."
        ))
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
        message(paste(
            "Variance-stabilizing transformation (vst) was not calculated.",
            "Skipping."
        ))
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
plotMeanSD.bcbioRNASeq <-  # nolint
    function(object, legend) {
        .plotMeanSD(
            raw = counts(object, normalized = FALSE),
            normalized = counts(object, normalized = TRUE),
            rlog = assays(object)[["rlog"]],
            vst = assays(object)[["vst"]],
            legend = legend
        )
    }

formals(plotMeanSD.bcbioRNASeq)[["legend"]] <- formalsList[["legend"]]



#' @rdname plotMeanSD
#' @export
setMethod(
    f = "plotMeanSD",
    signature = signature("bcbioRNASeq"),
    definition = plotMeanSD.bcbioRNASeq
)
