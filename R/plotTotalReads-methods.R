#' @name plotTotalReads
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inherit bioverbs::plotTotalReads
#' @note Updated 2019-09-15.
#'
#' @inheritParams acidroxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' data(bcb)
#' plotTotalReads(bcb)
NULL



#' @rdname plotTotalReads
#' @name plotTotalReads
#' @importFrom bioverbs plotTotalReads
#' @usage plotTotalReads(object, ...)
#' @export
NULL



## Updated 2019-09-15.
`plotTotalReads,bcbioRNASeq` <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        limit = 20e6L,
        perMillion = TRUE,
        fill,
        flip,
        title = "Total reads"
    ) {
        validObject(object)
        assert(
            isInt(limit),
            isNonNegative(limit),
            isFlag(perMillion),
            isGGScale(fill, scale = "discrete", aes = "fill", nullOK = TRUE),
            isFlag(flip),
            isString(title, nullOK = TRUE)
        )
        interestingGroups(object) <-
            matchInterestingGroups(object, interestingGroups)
        interestingGroups <- interestingGroups(object)
        data <- metrics(object)
        ## Convert to per million, if desired.
        yLab <- "reads"
        if (isTRUE(perMillion)) {
            data[["totalReads"]] <- data[["totalReads"]] / 1e6L
            yLab <- paste(yLab, "per million")
        }
        p <- ggplot(
                data = data,
                mapping = aes(
                    x = !!sym("sampleName"),
                    y = !!sym("totalReads"),
                    fill = !!sym("interestingGroups")
                )
            ) +
            acid_geom_bar() +
            acid_scale_y_continuous_nopad() +
            labs(
                title = title,
                x = NULL,
                y = yLab,
                fill = paste(interestingGroups, collapse = ":\n")
            )
        ## Limit.
        if (isPositive(limit)) {
            if (isTRUE(perMillion)) {
                if (limit < 1e6L) {
                    warning("'limit': Use absolute value, not per million.")
                } else {
                    limit <- limit / 1e6L
                }
            }
            p <- p + acid_geom_abline(yintercept = limit)
        }
        ## Fill.
        if (is(fill, "ScaleDiscrete")) {
            p <- p + fill
        }
        ## Flip.
        if (isTRUE(flip)) {
            p <- acid_coord_flip(p)
        }
        ## Hide sample name legend.
        if (identical(interestingGroups, "sampleName")) {
            p <- p + guides(fill = FALSE)
        }
        ## Return.
        p
    }

formals(`plotTotalReads,bcbioRNASeq`)[c("fill", "flip")] <-
    formalsList[c("fill.discrete", "flip")]



#' @rdname plotTotalReads
#' @export
setMethod(
    f = "plotTotalReads",
    signature = signature("bcbioRNASeq"),
    definition = `plotTotalReads,bcbioRNASeq`
)
