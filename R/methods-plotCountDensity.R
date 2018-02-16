#' Plot Count Density
#'
#' @rdname plotCountDensity
#' @name plotCountDensity
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inherit plotTotalReads
#' @inheritParams plotGene
#'
#' @param normalized Count normalization method. See [counts()] documentation
#'   for more information.
#' @param style Desired plot style (`line` or `solid`).
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq
#' plotCountDensity(bcb, style = "solid")
#' plotCountDensity(
#'     bcb,
#'     style = "line",
#'     interestingGroups = "sampleName",
#'     fill = NULL)
#'
#' # data.frame
#' df <- meltLog10(bcb, normalized = "tmm")
#' plotCountDensity(df)
NULL



# Constructors =================================================================
#' @importFrom bcbioBase uniteInterestingGroups
#' @importFrom ggplot2 aes_string geom_density ggplot guides labs
.plotCountDensity <- function(
    object,
    interestingGroups = "sampleName",
    style = "solid",
    color = scale_color_viridis(discrete = TRUE),
    fill = scale_fill_viridis(discrete = TRUE),
    title = TRUE) {
    assert_is_data.frame(object)
    assert_formal_interesting_groups(object, interestingGroups)
    assert_is_a_string(style)
    assert_is_subset(style, c("line", "solid"))
    assert_is_any_of(color, c("ScaleDiscrete", NULL))
    assert_is_any_of(fill, c("ScaleDiscrete", NULL))
    assert_is_any_of(title, c("character", "logical", "NULL"))

    data <- uniteInterestingGroups(object, interestingGroups)

    if (isTRUE(title)) {
        title <- "count density"
    } else if (!is.character(title)) {
        title <- NULL
    }

    p <- ggplot(
        data = data,
        mapping = aes_string(
            x = "counts",
            group = "interestingGroups",
            color = "interestingGroups",
            fill = "interestingGroups")
    ) +
        labs(
            title = title,
            x = "log10 counts per gene",
            fill = paste(interestingGroups, collapse = ":\n"))

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



.plotCountDensity.bcbioRNASeq <- function(  # nolint
    object,
    interestingGroups,
    normalized = "tmm",
    style = "solid",
    color = scale_color_viridis(discrete = TRUE),
    fill = scale_fill_viridis(discrete = TRUE),
    title = TRUE) {
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }
    .plotCountDensity(
        meltLog10(object, normalized = normalized),
        interestingGroups = interestingGroups,
        style = style,
        color = color,
        fill = fill,
        title = title)
}



# Methods ======================================================================
#' @rdname plotCountDensity
#' @export
setMethod(
    "plotCountDensity",
    signature("bcbioRNASeq"),
    .plotCountDensity.bcbioRNASeq)



#' @rdname plotCountDensity
#' @export
setMethod(
    "plotCountDensity",
    signature("data.frame"),
    .plotCountDensity)
