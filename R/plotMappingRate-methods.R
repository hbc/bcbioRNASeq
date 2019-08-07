#' @name plotMappingRate
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inherit bioverbs::plotMappingRate
#' @note Updated 2019-08-07.
#'
#' @inheritParams acidroxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' data(bcb)
#' plotMappingRate(bcb)
NULL



#' @rdname plotMappingRate
#' @name plotMappingRate
#' @importFrom bioverbs plotMappingRate
#' @usage plotMappingRate(object, ...)
#' @export
NULL



## Updated 2019-07-23.
`plotMappingRate,bcbioRNASeq` <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        limit = 0.7,
        fill,
        flip,
        title = "Mapping rate"
    ) {
        validObject(object)
        assert(
            isNumber(limit),
            isProportion(limit),
            isGGScale(fill, scale = "discrete", aes = "fill", nullOK = TRUE),
            isFlag(flip),
            isString(title, nullOK = TRUE)
        )
        interestingGroups(object) <-
            matchInterestingGroups(object, interestingGroups)
        interestingGroups <- interestingGroups(object)

        p <- ggplot(
            data = metrics(object),
            mapping = aes(
                x = !!sym("sampleName"),
                y = !!sym("mappedReads") / !!sym("totalReads") * 100L,
                fill = !!sym("interestingGroups")
            )
        ) +
            acid_geom_bar() +
            acid_scale_y_continuous_nopad(limits = c(0L, 100L)) +
            labs(
                title = title,
                x = NULL,
                y = "mapping rate (%)",
                fill = paste(interestingGroups, collapse = ":\n")
            )

        if (isPositive(limit)) {
            ## Convert to percentage.
            if (limit > 1L) {
                ## nocov start
                warning("`limit`: Use ratio (0-1) instead of percentage.")
                ## nocov end
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
            p <- acid_coord_flip(p)
        }

        if (identical(interestingGroups, "sampleName")) {
            p <- p + guides(fill = FALSE)
        }

        p
    }

formals(`plotMappingRate,bcbioRNASeq`)[c("fill", "flip")] <-
    formalsList[c("fill.discrete", "flip")]



#' @rdname plotMappingRate
#' @export
setMethod(
    f = "plotMappingRate",
    signature = signature("bcbioRNASeq"),
    definition = `plotMappingRate,bcbioRNASeq`
)
