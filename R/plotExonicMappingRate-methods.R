#' @name plotExonicMappingRate
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inherit bioverbs::plotExonicMappingRate
#' @inheritParams minimalism::params
#' @inheritParams basejump::params
#' @inheritParams params
#'
#' @description
#' Ideally, at least 60 percent of total reads should map to exons for RNA-seq.
#'
#' @examples
#' data(bcb)
#' plotExonicMappingRate(bcb)
NULL



#' @importFrom bioverbs plotExonicMappingRate
#' @aliases NULL
#' @export
bioverbs::plotExonicMappingRate



plotExonicMappingRate.bcbioRNASeq <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        limit = 0.6,
        fill,
        flip,
        title = "exonic mapping rate"
    ) {
        validObject(object)
        interestingGroups(object) <-
            matchInterestingGroups(object, interestingGroups)
        interestingGroups <- interestingGroups(object)
        assert(
            isNumber(limit),
            isProportion(limit),
            isGGScale(fill, scale = "discrete", aes = "fill", nullOK = TRUE),
            isFlag(flip),
            isString(title, nullOK = TRUE)
        )

        p <- metrics(object) %>%
            ggplot(
                mapping = aes(
                    x = !!sym("sampleName"),
                    y = !!sym("exonicRate") * 100L,
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
                y = "exonic mapping rate (%)",
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
                p <- p + acid_geom_abline(yintercept = limit)
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

formals(plotExonicMappingRate.bcbioRNASeq)[["fill"]] <-
    formalsList[["fill.discrete"]]
formals(plotExonicMappingRate.bcbioRNASeq)[["flip"]] <-
    formalsList[["flip"]]



#' @rdname plotExonicMappingRate
#' @export
setMethod(
    f = "plotExonicMappingRate",
    signature = signature("bcbioRNASeq"),
    definition = plotExonicMappingRate.bcbioRNASeq
)
