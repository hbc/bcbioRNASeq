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
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq
#' plotCountsPerGene(bcb)
#' plotCountsPerGene(
#'     bcb,
#'     interestingGroups = "sampleName",
#'     fill = NULL)
#'
#' # data.frame
#' df <- meltLog10(bcb, normalized = "tmm")
#' plotCountsPerGene(df)
NULL



# Constructors =================================================================
#' @importFrom bcbioBase uniteInterestingGroups
#' @importFrom ggplot2 aes_string geom_boxplot ggplot guides labs
.plotCountsPerGene <- function(
    object,
    interestingGroups = "sampleName",
    fill = scale_fill_viridis(discrete = TRUE),
    flip = TRUE,
    title = TRUE) {
    assert_is_data.frame(object)
    assert_formal_interesting_groups(object, interestingGroups)
    .assert_formal_discrete_scale(fill)
    assert_is_a_bool(flip)
    .assert_formal_title(title)
    
    data <- uniteInterestingGroups(object, interestingGroups)

    if (isTRUE(title)) {
        title <- "counts per gene"
    } else if (!is.character(title)) {
        title <- NULL
    }

    p <- ggplot(
        data = data,
        mapping = aes_string(
            x = "sampleName",
            y = "counts",
            fill = "interestingGroups")
    ) +
        geom_boxplot(color = lineColor, outlier.shape = NA) +
        labs(
            title = title,
            x = "sample",
            y = "log10 counts per gene",
            fill = paste(interestingGroups, collapse = ":\n"))

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



.plotCountsPerGene.bcbioRNASeq <- function(  # nolint
    object,
    interestingGroups,
    normalized = "tmm",
    fill = scale_fill_viridis(discrete = TRUE),
    flip = TRUE,
    title = TRUE) {
    # Passthrough: fill, flip, title
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }
    assert_is_a_string(normalized)
    .plotCountsPerGene(
        meltLog10(object, normalized = normalized),
        interestingGroups = interestingGroups,
        fill = fill,
        flip = flip,
        title = title)
}



# Methods ======================================================================
#' @rdname plotCountsPerGene
#' @export
setMethod(
    "plotCountsPerGene",
    signature("bcbioRNASeq"),
    .plotCountsPerGene.bcbioRNASeq)



#' @rdname plotCountsPerGene
#' @export
setMethod(
    "plotCountsPerGene",
    signature("data.frame"),
    .plotCountsPerGene)
