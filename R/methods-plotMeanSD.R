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
#' @param orientation Orientation to use for plot grid, either `horizontal` or
#'   `vertical`.
#' @param legend Include the color bar legend. This is typically not that
#'   informative and is disabled by default, to improve the plot appearance.
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
    legend = FALSE,
    headerLevel = 2L,
    return = c("grid", "list", "markdown")
) {
    assert_is_matrix(raw)
    assert_is_matrix(normalized)
    assert_is_matrix(rlog)
    assert_is_matrix(vst)
    assert_is_a_bool(legend)
    return <- match.arg(return)

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

    # Remove the plot (color) legend, if desired
    if (!isTRUE(legend)) {
        gglog2 <- gglog2 +
            theme(legend.position = "none")
        ggrlog <- ggrlog +
            theme(legend.position = "none")
        ggvst <- ggvst +
            theme(legend.position = "none")
    }

    plotlist <- list(
        log2 = gglog2,
        rlog = ggrlog,
        vst = ggvst
    )

    dynamicPlotlist(
        plotlist = plotlist,
        headerLevel = headerLevel,
        return = return
    )
}


# Methods ======================================================================
#' @rdname plotMeanSD
#' @export
setMethod(
    "plotMeanSD",
    signature("bcbioRNASeq"),
    function(
        object,
        legend = FALSE,
        headerLevel = 2L,
        return = c("grid", "list", "markdown")
    ) {
        # Passthrough: legend, headerLevel
        return <- match.arg(return)

        # Require that the DESeq2 transformations are slotted.
        # If `transformationLimit` was applied, this function will error.
        rlog <- assays(object)[["rlog"]]
        vst <- assays(object)[["vst"]]
        assert_is_matrix(rlog)
        assert_is_matrix(vst)

        .plotMeanSD(
            raw = counts(object, normalized = FALSE),
            normalized = counts(object, normalized = TRUE),
            rlog = rlog,
            vst = vst,
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
        legend = FALSE,
        headerLevel = 2L,
        return = c("grid", "list", "markdown")
    ) {
        # Passthrough: legend, headerLevel
        return <- match.arg(return)

        raw <- counts(object, normalized = FALSE)
        normalized <- counts(object, normalized = TRUE)
        rlog <- assay(rlog(object))
        vst <- assay(varianceStabilizingTransformation(object))

        .plotMeanSD(
            raw = raw,
            normalized = normalized,
            rlog = rlog,
            vst = vst,
            legend = legend,
            headerLevel = headerLevel,
            return = return
        )
    }
)
