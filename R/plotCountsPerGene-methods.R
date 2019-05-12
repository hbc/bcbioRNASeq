#' @name plotCountsPerGene
#' @inherit bioverbs::plotCountsPerGene
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @section TMM:
#' We recommend visualizing counts normalized with the Trimmed Mean of M-Values
#' (TMM) method here (Robinson, et al., 2010). TMM normalization equates the
#' overall expression levels of genes between samples under the assumption that
#' the majority of them are not differentially expressed. Therefore, by
#' normalizing for total RNA expression by sample, we expect the spread of the
#' TMM-normalized counts per gene to be similar for every sample.
#'
#' @inheritParams general
#' @param ... Additional arguments.
#'
#' @return `ggplot`.
#'
#' @examples
#' plotCountsPerGene(bcb_small)
NULL



#' @rdname plotCountsPerGene
#' @name plotCountsPerGene
#' @importFrom bioverbs plotCountsPerGene
#' @usage plotCountsPerGene(object, ...)
#' @export
NULL



plotCountsPerGene.bcbioRNASeq <-  # nolint
    function(
        object,
        interestingGroups,
        normalized = c("tmm", "vst", "rlog", "tpm", "rle"),
        fill = getOption("bcbio.discrete.fill", NULL),
        flip = getOption("bcbio.flip", TRUE),
        title = "counts per gene"
    ) {
        # Passthrough: fill, flip, title
        validObject(object)
        interestingGroups <- matchInterestingGroups(
            object = object,
            interestingGroups = interestingGroups
        )
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
            mapping = aes(
                x = !!sym("sampleName"),
                y = !!sym("counts"),
                fill = !!sym("interestingGroups")
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



#' @rdname plotCountsPerGene
#' @export
setMethod(
    f = "plotCountsPerGene",
    signature = signature("bcbioRNASeq"),
    definition = plotCountsPerGene.bcbioRNASeq
)
