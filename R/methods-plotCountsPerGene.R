#' Plot Counts Per Gene
#'
#' @name plotCountsPerGene
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' plotCountsPerGene(bcb_small)
NULL



# Methods ======================================================================
#' @rdname plotCountsPerGene
#' @export
setMethod(
    "plotCountsPerGene",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        normalized = c("tmm", "vst", "rlog", "tpm", "rle"),
        fill = NULL,
        flip = TRUE,
        title = "counts per gene"
    ) {
        # Passthrough: fill, flip, title
        validObject(object)
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        } else {
            interestingGroups(object) <- interestingGroups
        }
        normalized <- match.arg(normalized)
        assertIsFillScaleDiscreteOrNULL(fill)
        assert_is_a_bool(flip)
        assertIsAStringOrNULL(title)

        # Subset the counts matrix to only include non-zero genes
        nonzero <- .nonzeroGenes(object)
        counts <- counts(object, normalized = normalized)
        counts <- counts[nonzero, , drop = FALSE]

        # Apply log2 transformation, if  necessary
        if (normalized %in% c("rlog", "vst")) {
            # Already log2
            fxn <- .meltCounts
        } else {
            fxn <- .meltLog2Counts
        }

        data <- fxn(counts, sampleData = sampleData(object))

        # Subtitle
        if (is_a_string(title)) {
            subtitle <- paste(nrow(counts), "non-zero genes")
        } else {
            subtitle <- NULL
        }

        p <- ggplot(
            data = data,
            mapping = aes_string(
                x = "sampleName",
                y = "counts",
                fill = "interestingGroups"
            )
        ) +
            geom_boxplot(color = "black", outlier.shape = NA) +
            labs(
                title = title,
                subtitle = subtitle,
                x = NULL,
                y = paste(normalized, "counts (log2)"),
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
)
