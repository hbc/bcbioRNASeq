#' Plot Counts Per Gene
#'
#' @rdname plotCountsPerGene
#' @name plotCountsPerGene
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inherit plotTotalReads
#' @inheritParams plotCountDensity
#'
#' @examples
#' plotCountsPerGene(bcb_small)
NULL



# Constructors =================================================================
#' @importFrom bcbioBase uniteInterestingGroups
#' @importFrom ggplot2 aes_string geom_boxplot ggplot guides labs
.plotCountsPerGene <- function(
    object,
    interestingGroups = "sampleName",
    fill = scale_fill_viridis(discrete = TRUE),
    flip = TRUE,
    title = TRUE,
    xlab = "counts"
) {
    assert_is_data.frame(object)
    assertFormalInterestingGroups(object, interestingGroups)
    assertIsFillScaleDiscreteOrNULL(fill)
    assert_is_a_bool(flip)

    # Title
    if (isTRUE(title)) {
        title <- "counts per gene"
    } else if (!is_a_string(title)) {
        title <- NULL
    }

    data <- uniteInterestingGroups(object, interestingGroups)

    p <- ggplot(
        data = data,
        mapping = aes_string(
            x = "sampleName",
            y = "counts",
            fill = "interestingGroups"
        )
    ) +
        geom_boxplot(color = lineColor, outlier.shape = NA) +
        labs(
            title = title,
            x = "sample",
            y = xlab,
            fill = paste(interestingGroups, collapse = ":\n")
        )

    if (is(fill, "ScaleDiscrete")) {
        p <- p + fill
    }

    if (isTRUE(flip)) {
        p <- p + coord_flip()
    }

    if (identical(interestingGroups, "sampleName")) {
        p <- p + guides(fill = FALSE)
    }

    p
}



# Methods ======================================================================
#' @rdname plotCountsPerGene
#' @export
setMethod(
    "plotCountsPerGene",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        normalized = "rlog",
        fill = scale_fill_viridis(discrete = TRUE),
        flip = TRUE,
        title = TRUE
    ) {
        # Passthrough: fill, flip, title
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        assert_is_a_string(normalized)
        if (normalized %in% c("rlog", "vst")) {
            fxn <- .meltCounts
        } else {
            fxn <- .meltLog2Counts
        }
        melt <- fxn(
            counts(object, normalized = normalized),
            colData = colData(object)
        )
        .plotCountsPerGene(
            object = melt,
            interestingGroups = interestingGroups,
            fill = fill,
            flip = flip,
            title = title,
            xlab = paste("log2", normalized, "counts")
        )
    })
