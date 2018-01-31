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
#' @importFrom viridis scale_fill_viridis
.plotRRNAMappingRate <- function(
    object,
    interestingGroups = "sampleName",
    warnLimit = 10L,
    fill = viridis::scale_fill_viridis(discrete = TRUE),
    flip = TRUE,
    title = TRUE) {
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

    metrics <- uniteInterestingGroups(object, interestingGroups)
    p <- ggplot(
        metrics,
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
#' @rdname plotRRNAMappingRate
#' @importFrom viridis scale_fill_viridis
#' @export
setMethod(
    "plotRRNAMappingRate",
    signature("bcbioRNASeq"),
    function(
        object,
        interestingGroups,
        warnLimit = 10L,
        fill = viridis::scale_fill_viridis(discrete = TRUE),
        flip = TRUE,
        title = TRUE) {
        if (is.null(metrics(object))) return(NULL)
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        .plotRRNAMappingRate(
            metrics(object),
            interestingGroups = interestingGroups,
            warnLimit = warnLimit,
            fill = fill,
            flip = flip,
            title = title)
    })



#' @rdname plotRRNAMappingRate
#' @export
setMethod(
    "plotRRNAMappingRate",
    signature("data.frame"),
    .plotRRNAMappingRate)
