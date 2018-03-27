#' Plot Count Density
#'
#' @name plotCountDensity
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inheritParams general
#' @param style Desired plot style ("`line`" or "`solid`").
#'
#' @return `ggplot`.
#'
#' @examples
#' # bcbioRNASeq ====
#' plotCountDensity(bcb_small)
#' plotCountDensity(
#'     object = bcb_small,
#'     style = "line",
#'     interestingGroups = "sampleName"
#' )
NULL



# Constructors =================================================================
.plotCountDensity <- function(
    object,
    interestingGroups,
    normalized = c("tmm", "rlog", "vst", "tpm"),
    style = c("line", "solid"),
    color = scale_color_viridis(discrete = TRUE),
    fill = scale_fill_viridis(discrete = TRUE),
    title = "count density"
) {
    validObject(object)
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }
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
    data <- fxn(counts, colData = colData(object)) %>%
        uniteInterestingGroups(interestingGroups)

    # Subtitle
    if (is_a_string(title)) {
        subtitle <- paste(nrow(counts), "non-zero genes")
    } else {
        subtitle <- NULL
    }

    p <- ggplot(
        data = data,
        mapping = aes_string(
            x = "counts",
            group = "interestingGroups",
            color = "interestingGroups",
            fill = "interestingGroups"
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

    if (identical(interestingGroups, "sampleName")) {
        p <- p + guides(fill = FALSE)
    }

    p
}



# Methods ======================================================================
#' @rdname plotCountDensity
#' @export
setMethod(
    "plotCountDensity",
    signature("bcbioRNASeq"),
    .plotCountDensity
)
