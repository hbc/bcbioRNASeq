#' @name plotMappingRate
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inherit basejump::plotMappingRate
#' @inheritParams basejump::params
#' @inheritParams params
#' @examples
#' data(bcb)
#' plotMappingRate(bcb)
NULL



#' @importFrom basejump plotMappingRate
#' @aliases NULL
#' @export
basejump::plotMappingRate



plotMappingRate.bcbioRNASeq <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        limit = 0.7,
        fill,
        flip,
        title = "mapping rate"
    ) {
        validObject(object)
        interestingGroups(object) <-
            matchInterestingGroups(object, interestingGroups)
        assert(
            isNumber(limit),
            isProportion(limit),
            isGGScale(fill, scale = "discrete", aes = "fill") || is.null(fill),
            isFlag(flip),
            isString(title) || is.null(title)
        )

        p <- ggplot(
            data = metrics(object),
            mapping = aes(
                x = !!sym("sampleName"),
                y = !!sym("mappedReads") / !!sym("totalReads") * 100L,
                fill = !!sym("interestingGroups")
            )
        ) +
            geom_bar(
                color = "black",
                stat = "identity"
            ) +
            ylim(0L, 100L) +
            labs(
                title = title,
                x = NULL,
                y = "mapping rate (%)",
                fill = paste(interestingGroups, collapse = ":\n")
            )

        if (isPositive(limit)) {
            # Convert to percentage.
            if (limit > 1L) {
                # nocov start
                warning("`limit`: Use ratio (0-1) instead of percentage.")
                # nocov end
            } else {
                limit <- limit * 100L
            }
            if (limit < 100L) {
                p <- p + basejump_geom_abline(yintercept = limit)
            }
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

formals(plotMappingRate.bcbioRNASeq)[["fill"]] <-
    formalsList[["fill.discrete"]]
formals(plotMappingRate.bcbioRNASeq)[["flip"]] <-
    formalsList[["flip"]]



#' @rdname plotMappingRate
#' @export
setMethod(
    f = "plotMappingRate",
    signature = signature("bcbioRNASeq"),
    definition = plotMappingRate.bcbioRNASeq
)
