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
.plotMeanSD.assays <- function(  # nolint
    raw,
    normalized,
    rlog,
    vst,
    orientation = c("vertical", "horizontal"),
    legend = FALSE
) {
    assert_is_matrix(raw)
    assert_is_matrix(normalized)
    assert_is_matrix(rlog)
    assert_is_matrix(vst)
    orientation <- match.arg(orientation)
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

    # Return either horizontal or vertical
    if (orientation == "horizontal") {
        ncol <- 3L
        nrow <- 1L
    } else if (orientation == "vertical") {
        ncol <- 1L
        nrow <- 3L
    }

    plot_grid(
        gglog2,
        ggrlog,
        ggvst,
        labels = "AUTO",
        ncol = ncol,
        nrow = nrow
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
        orientation = c("vertical", "horizontal"),
        legend = FALSE
    ) {
        # Passthrough: orientation, legend
        # Require that the DESeq2 transformations are slotted
        rlog <- assays(object)[["rlog"]]
        vst <- assays(object)[["vst"]]
        assert_is_matrix(rlog)
        assert_is_matrix(vst)
        .plotMeanSD.assays(
            raw = counts(object, normalized = FALSE),
            normalized = counts(object, normalized = TRUE),
            rlog = rlog,
            vst = vst,
            orientation = orientation,
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
        orientation = c("vertical", "horizontal"),
        legend = FALSE
    ) {
        # Passthrough arguments: orientation, legend
        .plotMeanSD.assays(
            raw = counts(object, normalized = FALSE),
            normalized = counts(object, normalized = TRUE),
            rlog = assay(rlog(object)),
            vst = assay(varianceStabilizingTransformation(object)),
            orientation = orientation,
            legend = legend
        )
    }
)
