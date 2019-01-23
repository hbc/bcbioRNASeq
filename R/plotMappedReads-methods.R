#' @name plotMappedReads
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inherit bioverbs::plotMappedReads
#' @inheritParams basejump::params
#' @inheritParams params
#' @examples
#' data(bcb)
#' plotMappedReads(bcb)
NULL



#' @importFrom bioverbs plotMappedReads
#' @aliases NULL
#' @export
bioverbs::plotMappedReads



plotMappedReads.bcbioRNASeq <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        limit = 10e6L,
        fill,
        flip,
        title = "mapped reads"
    ) {
        validObject(object)
        interestingGroups(object) <-
            matchInterestingGroups(object, interestingGroups)
        interestingGroups <- interestingGroups(object)
        assert(
            isInt(limit),
            isNonNegative(limit),
            isGGScale(fill, scale = "discrete", aes = "fill", nullOK = TRUE),
            isFlag(flip),
            isString(title, nullOK = TRUE)
        )

        p <- ggplot(
            data = metrics(object),
            mapping = aes(
                x = !!sym("sampleName"),
                y = !!sym("mappedReads") / 1e6L,
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
                y = "mapped reads per million",
                fill = paste(interestingGroups, collapse = ":\n")
            )

        if (isPositive(limit)) {
            # Convert limit to per million.
            if (limit < 1e6L) {
                # nocov start
                warning("`limit`: Use absolute value, not per million.")
                # nocov end
            } else {
                limit <- limit / 1e6L
            }
            if (limit >= 1L) {
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

formals(plotMappedReads.bcbioRNASeq)[["fill"]] <-
    formalsList[["fill.discrete"]]
formals(plotMappedReads.bcbioRNASeq)[["flip"]] <-
    formalsList[["flip"]]



#' @rdname plotMappedReads
#' @export
setMethod(
    f = "plotMappedReads",
    signature = signature("bcbioRNASeq"),
    definition = plotMappedReads.bcbioRNASeq
)
