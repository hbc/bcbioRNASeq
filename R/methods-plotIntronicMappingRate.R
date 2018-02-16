#' Plot Intronic Mapping Rate
#'
#' @rdname plotIntronicMappingRate
#' @name plotIntronicMappingRate
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
#' plotIntronicMappingRate(bcb)
#' plotIntronicMappingRate(
#'     bcb,
#'     interestingGroups = "sampleName",
#'     fill = NULL)
#'
#' # data.frame
#' df <- metrics(bcb)
#' plotIntronicMappingRate(df)
NULL



# Constructors =================================================================
#' @importFrom bcbioBase uniteInterestingGroups
#' @importFrom ggplot2 aes_ coord_flip geom_bar ggplot guides labs ylim
.plotIntronicMappingRate <- function(
    object,
    interestingGroups = "sampleName",
    warnLimit = 20L,
    fill = scale_fill_viridis(discrete = TRUE),
    flip = TRUE,
    title = TRUE) {
    if (isTRUE(title)) {
        title <- "intronic mapping rate"
    } else if (!is.character(title)) {
        title <- NULL
    }

    metrics <- uniteInterestingGroups(object, interestingGroups)
    p <- ggplot(
        metrics,
        mapping = aes_(
            x = ~sampleName,
            y = ~intronicRate * 100L,
            fill = ~interestingGroups)
    ) +
        geom_bar(stat = "identity") +
        labs(
            title = title,
            x = "sample",
            y = "intronic mapping rate (%)",
            fill = paste(interestingGroups, collapse = ":\n")
        ) +
        ylim(0L, 100L)

    if (is.numeric(warnLimit)) {
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



# Methods ======================================================================
#' @rdname plotIntronicMappingRate
#' @export
setMethod(
    "plotIntronicMappingRate",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        warnLimit = 20L,
        fill = scale_fill_viridis(discrete = TRUE),
        flip = TRUE,
        title = TRUE) {
        if (is.null(metrics(object))) return(NULL)
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        .plotIntronicMappingRate(
            metrics(object),
            interestingGroups = interestingGroups,
            warnLimit = warnLimit,
            fill = fill,
            flip = flip,
            title = title)
    })



#' @rdname plotIntronicMappingRate
#' @export
setMethod(
    "plotIntronicMappingRate",
    signature("data.frame"),
    .plotIntronicMappingRate)
