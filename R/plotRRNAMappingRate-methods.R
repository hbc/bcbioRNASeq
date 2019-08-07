#' @name plotRRNAMappingRate
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inherit bioverbs::plotRRNAMappingRate
#'
#' @inheritParams acidroxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' data(bcb)
#' plotRRNAMappingRate(bcb)
NULL



#' @rdname plotRRNAMappingRate
#' @name plotRRNAMappingRate
#' @importFrom bioverbs plotRRNAMappingRate
#' @usage plotRRNAMappingRate(object, ...)
#' @export
NULL



## Updated 2019-07-23.
`plotRRNAMappingRate,bcbioRNASeq` <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        limit = 0.1,
        fill,
        flip,
        title = "rRNA mapping rate"
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

        data <- metrics(object)

        ## Warn and early return if rRNA rate was not calculated.
        if (!"rrnaRate" %in% colnames(data)) {
            warning("rRNA mapping rate was not calculated. Skipping plot.")
            return(invisible())
        }

        p <- ggplot(
            data = data,
            mapping = aes(
                x = !!sym("sampleName"),
                y = !!sym("rrnaRate") * 100L,
                fill = !!sym("interestingGroups")
            )
        ) +
            acid_geom_bar() +
            acid_scale_y_continuous_nopad(limits = c(0L, 100L)) +
            labs(
                title = title,
                x = NULL,
                y = "rRNA mapping rate (%)",
                fill = paste(interestingGroups, collapse = ":\n")
            )

        if (isPositive(limit)) {
            ## Convert to percentage
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

formals(`plotRRNAMappingRate,bcbioRNASeq`)[c("fill", "flip")] <-
    formalsList[c("fill.discrete", "flip")]



#' @rdname plotRRNAMappingRate
#' @export
setMethod(
    f = "plotRRNAMappingRate",
    signature = signature("bcbioRNASeq"),
    definition = `plotRRNAMappingRate,bcbioRNASeq`
)
