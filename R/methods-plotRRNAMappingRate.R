#' Plot Ribosomal RNA (rRNA) Mapping Rate
#'
#' @rdname plotRRNAMappingRate
#' @name plotRRNAMappingRate
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inherit plotTotalReads
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioRNASeq"))
#'
#' # bcbioRNASeq
#' plotRRNAMappingRate(bcb)
#' plotRRNAMappingRate(
#'     bcb,
#'     interestingGroups = "sampleName",
#'     fill = NULL,
#'     warnLimit = NULL)
#'
#' # data.frame
#' df <- metrics(bcb)
#' plotRRNAMappingRate(df)
NULL



# Constructors =================================================================
#' @importFrom bcbioBase uniteInterestingGroups
#' @importFrom ggplot2 aes_ coord_flip geom_bar ggplot guides labs
.plotRRNAMappingRate <- function(
    object,
    interestingGroups = "sampleName",
    warnLimit = 10L,
    fill = scale_fill_viridis(discrete = TRUE),
    flip = TRUE,
    title = TRUE) {
    assert_is_data.frame(object)
    assert_formal_interesting_groups(object, interestingGroups)
    assert_is_an_implicit_integer(warnLimit)
    assert_all_are_non_negative(warnLimit)
    .assert_formal_scale_discrete(fill)
    assert_is_a_bool(flip)
    .assert_formal_title(title)

    if (isTRUE(title)) {
        title <- "rRNA mapping rate"
    } else if (!is.character(title)) {
        title <- NULL
    }

    # Fix for legacy camel variant mismatch (e.g. rRnaRate).
    if (!"rrnaRate" %in% colnames(object)) {
        # grep match the outdated camel variant
        col <- grep(
            x = colnames(object),
            pattern = "rrnarate",
            ignore.case = TRUE,
            value = TRUE)
        object[["rrnaRate"]] <- object[[col]]
        object[[col]] <- NULL
    }

    data <- uniteInterestingGroups(object, interestingGroups)

    p <- ggplot(
        data = data,
        mapping = aes_(
            x = ~sampleName,
            y = ~rrnaRate * 100L,
            fill = ~interestingGroups)
    ) +
        geom_bar(stat = "identity") +
        labs(
            title = title,
            x = "sample",
            y = "rRNA mapping rate (%)",
            fill = paste(interestingGroups, collapse = ":\n")
        )

    if (is_positive(warnLimit)) {
        p <- p + qcWarnLine(warnLimit)
    }

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



.plotRRNAMappingRate.bcbioRNASeq <- function(  # nolint
    object,
    interestingGroups,
    warnLimit = 10L,
    fill = scale_fill_viridis(discrete = TRUE),
    flip = TRUE,
    title = TRUE) {
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }
    .plotRRNAMappingRate(
        object = metrics(object),
        interestingGroups = interestingGroups,
        warnLimit = warnLimit,
        fill = fill,
        flip = flip,
        title = title)
}



# Methods ======================================================================
#' @rdname plotRRNAMappingRate
#' @export
setMethod(
    "plotRRNAMappingRate",
    signature("bcbioRNASeq"),
    .plotRRNAMappingRate.bcbioRNASeq)



#' @rdname plotRRNAMappingRate
#' @export
setMethod(
    "plotRRNAMappingRate",
    signature("data.frame"),
    .plotRRNAMappingRate)
