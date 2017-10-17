#' Plot Ribosomal RNA (rRNA) Mapping Rate
#'
#' @rdname plotRRNAMappingRate
#' @name plotRRNAMappingRate
#' @family Quality Control Plots
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inherit qcPlots
#'
#' @inheritParams AllGenerics
#'
#' @examples
#' data(bcb)
#'
#' # bcbioRNASeq
#' plotRRNAMappingRate(bcb)
#'
#' \dontrun{
#' plotRRNAMappingRate(bcb, interestingGroups = "group")
#'
#' # data.frame
#' metrics(bcb) %>%
#'     plotRRNAMappingRate()
#' }
NULL



# Constructors ====
.plotRRNAMappingRate <- function(
    object,
    interestingGroups = "sampleName",
    warnLimit = 10,
    flip = TRUE) {
    # Fix for camel variant mismatch (e.g. rRnaRate).
    if (!"rrnaRate" %in% colnames(object)) {
        # grep match the outdated camel variant
        col <- grep(x = colnames(object),
                    pattern = "rrnarate",
                    ignore.case = TRUE,
                    value = TRUE)
        object[["rrnaRate"]] <- object[[col]]
        object[[col]] <- NULL
    }

    p <- ggplot(
        object,
        mapping = aes_(
            x = ~sampleName,
            y = ~rrnaRate * 100,
            fill = as.name(interestingGroups))
    ) +
        geom_bar(stat = "identity") +
        labs(title = "rrna mapping rate",
             x = "sample",
             y = "rRNA mapping rate (%)") +
        scale_fill_viridis(discrete = TRUE)
    if (!is.null(warnLimit)) {
        p <- p +
            qcWarnLine(warnLimit)
    }
    if (isTRUE(flip)) {
        p <- p +
            coord_flip()
    }
    p
}



# Methods ====
#' @rdname plotRRNAMappingRate
#' @export
setMethod(
    "plotRRNAMappingRate",
    signature("bcbioRNASeqANY"),
    function(
        object,
        interestingGroups,
        warnLimit = 10,
        flip = TRUE) {
        if (is.null(metrics(object))) {
            return(NULL)
        }
        if (missing(interestingGroups)) {
            interestingGroups <-
                metadata(object)[["interestingGroups"]][[1]]
        }
        .plotRRNAMappingRate(
            metrics(object),
            interestingGroups = interestingGroups,
            warnLimit = warnLimit,
            flip = flip)
    })



#' @rdname plotRRNAMappingRate
#' @export
setMethod(
    "plotRRNAMappingRate",
    signature("data.frame"),
    .plotRRNAMappingRate)
