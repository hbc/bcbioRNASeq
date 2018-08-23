#' Plot Count Density
#'
#' Generally, we expect similar count spreads for all genes between samples
#' unless the total expressed RNA per sample is different.
#'
#' @name plotCountDensity
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inheritParams general
#' @param style `string`. Desired plot style ("`line`" or "`solid`").
#'
#' @return `ggplot`.
#'
#' @examples
#' # bcbioRNASeq ====
#' plotCountDensity(bcb_small)
NULL



#' @rdname plotCountDensity
#' @export
setMethod(
    "plotCountDensity",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        normalized = c("tmm", "vst", "rlog", "tpm", "rle"),
        style = c("line", "solid"),
        color = getOption("bcbio.discrete.color", NULL),
        fill = getOption("bcbio.discrete.fill", NULL),
        title = "count density"
    ) {
        validObject(object)
        interestingGroups <- matchInterestingGroups(
            object = object,
            interestingGroups = interestingGroups
        )
        normalized <- match.arg(normalized)
        style <- match.arg(style)
        assertIsColorScaleDiscreteOrNULL(color)
        assertIsFillScaleDiscreteOrNULL(fill)
        assertIsAStringOrNULL(title)

        styleLab <- paste(interestingGroups, collapse = ":\n")

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

        # Melt the counts into long format
        sampleData <- sampleData(object, interestingGroups = interestingGroups)
        data <- fxn(counts, sampleData = sampleData)

        # Subtitle
        if (is_a_string(title)) {
            subtitle <- paste(nrow(counts), "non-zero genes")
        } else {
            subtitle <- NULL
        }

        p <- ggplot(
            data = data,
            mapping = aes(
                x = !!sym("counts"),
                group = !!sym("interestingGroups"),
                color = !!sym("interestingGroups"),
                fill = !!sym("interestingGroups")
            )
        ) +
            labs(
                title = title,
                subtitle = subtitle,
                x = paste(normalized, "counts (log2)"),
                color = styleLab,
                fill = styleLab
            )

        if (style == "line") {
            p <- p + geom_density(fill = NA)
            if (is(color, "ScaleDiscrete")) {
                p <- p + color
            }
        } else if (style == "solid") {
            p <- p + geom_density(alpha = 0.75, color = NA)
            if (is(fill, "ScaleDiscrete")) {
                p <- p + fill
            }
        }

        p
    }
)
