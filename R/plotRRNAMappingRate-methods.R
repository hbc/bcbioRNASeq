#' @name plotRRNAMappingRate
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inherit bioverbs::plotRRNAMappingRate
#' @inheritParams basejump::params
#' @inheritParams params
#' @examples
#' data(bcb)
#' plotRRNAMappingRate(bcb)
NULL



#' @importFrom bioverbs plotRRNAMappingRate
#' @aliases NULL
#' @export
bioverbs::plotRRNAMappingRate



plotRRNAMappingRate.bcbioRNASeq <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        limit = 0.1,
        fill,
        flip,
        title = "rRNA mapping rate"
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

        data <- metrics(object)

        # Warn and early return if rRNA rate was not calculated.
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
            geom_bar(
                color = "black",
                stat = "identity"
            ) +
            labs(
                title = title,
                x = NULL,
                y = "rRNA mapping rate (%)",
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

formals(plotRRNAMappingRate.bcbioRNASeq)[["fill"]] <-
    formalsList[["fill.discrete"]]
formals(plotRRNAMappingRate.bcbioRNASeq)[["flip"]] <-
    formalsList[["flip"]]



#' @rdname plotRRNAMappingRate
#' @export
setMethod(
    f = "plotRRNAMappingRate",
    signature = signature("bcbioRNASeq"),
    definition = plotRRNAMappingRate.bcbioRNASeq
)
