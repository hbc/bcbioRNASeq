# FIXME Add the number of non-zero genes to plot

#' Plot Counts Per Gene
#'
#' @name plotCountsPerGene
#' @family Quality Control Functions
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
.plotCountsPerGene.melted <- function(
    object,
    interestingGroups = "sampleName",
    fill = scale_fill_viridis(discrete = TRUE),
    flip = TRUE,
    title = NULL,
    subtitle = NULL,
    xlab = "counts"
) {
    assert_is_data.frame(object)
    assertFormalInterestingGroups(object, interestingGroups)
    assertIsFillScaleDiscreteOrNULL(fill)
    assert_is_a_bool(flip)
    assertIsAStringOrNULL(title)
    if (is.null(title)) {
        subtitle <- NULL
    }
    assertIsAStringOrNULL(subtitle)
    assert_is_a_string(xlab)

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
            subtitle = subtitle,
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
        normalized = c("rlog", "vst", "tmm", "tpm"),
        fill = scale_fill_viridis(discrete = TRUE),
        flip = TRUE,
        title = "counts per gene"
    ) {
        # Passthrough: fill, flip, title
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        normalized <- match.arg(normalized)

        # Subset the counts matrix to only include non-zero genes
        raw <- assay(object)
        nonzero <- rowSums(raw) > 0L
        genes <- rownames(raw[nonzero, , drop = FALSE])
        counts <- counts(object, normalized = normalized)
        counts <- counts[genes, , drop = FALSE]

        # Apply log2 transformation, if  necessary
        if (normalized %in% c("rlog", "vst")) {
            # Already log2
            fxn <- .meltCounts
        } else {
            fxn <- .meltLog2Counts
        }
        melted <- fxn(counts, colData = colData(object))

        .plotCountsPerGene.melted(
            object = melted,
            interestingGroups = interestingGroups,
            fill = fill,
            flip = flip,
            title = title,
            subtitle = paste(nrow(counts), "non-zero genes"),
            xlab = paste("log2", normalized, "counts")
        )
    }
)
