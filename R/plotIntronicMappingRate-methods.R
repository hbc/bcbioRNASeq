#' @name plotIntronicMappingRate
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inherit basejump::plotIntronicMappingRate
#' @inheritParams basejump::params
#' @inheritParams params
#' @examples
#' data(bcb)
#' plotIntronicMappingRate(bcb)
NULL



#' @importFrom basejump plotIntronicMappingRate
#' @aliases NULL
#' @export
basejump::plotIntronicMappingRate



plotIntronicMappingRate.bcbioRNASeq <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        limit = 0.2,
        fill,
        flip,
        title = "intronic mapping rate"
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
                y = !!sym("intronicRate") * 100L,
                fill = !!sym("interestingGroups")
            )
        ) +
            geom_bar(
                color = "black",
                stat = "identity"
            ) +
            labs(
                title = title,
                x = NULL,
                y = "intronic mapping rate (%)",
                fill = paste(interestingGroups, collapse = ":\n")
            )

        if (isPositive(limit)) {
            # Convert to percentage
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

formals(plotIntronicMappingRate.bcbioRNASeq)[["fill"]] <-
    formalsList[["fill.discrete"]]
formals(plotIntronicMappingRate.bcbioRNASeq)[["flip"]] <-
    formalsList[["flip"]]



#' @rdname plotIntronicMappingRate
#' @export
setMethod(
    f = "plotIntronicMappingRate",
    signature = signature("bcbioRNASeq"),
    definition = plotIntronicMappingRate.bcbioRNASeq
)
