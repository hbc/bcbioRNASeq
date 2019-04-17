#' @name plotTotalReads
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inherit bioverbs::plotTotalReads
#' @inheritParams minimalism::params
#' @inheritParams basejump::params
#' @inheritParams params
#' @examples
#' data(bcb)
#' plotTotalReads(bcb)
NULL



#' @rdname plotTotalReads
#' @name plotTotalReads
#' @importFrom bioverbs plotTotalReads
#' @export
NULL



plotTotalReads.bcbioRNASeq <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        limit = 20e6L,
        perMillion = TRUE,
        fill,
        flip,
        title = "total reads"
    ) {
        validObject(object)
        interestingGroups(object) <-
            matchInterestingGroups(object, interestingGroups)
        interestingGroups <- interestingGroups(object)
        assert(
            isInt(limit),
            isNonNegative(limit),
            isFlag(perMillion),
            isGGScale(fill, scale = "discrete", aes = "fill", nullOK = TRUE),
            isFlag(flip),
            isString(title, nullOK = TRUE)
        )

        data <- metrics(object)

        # Convert to per million, if desired.
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
            geom_bar(
                color = "black",
                stat = "identity"
            ) +
            labs(
                title = title,
                x = NULL,
                y = yLab,
                fill = paste(interestingGroups, collapse = ":\n")
            )

        if (isPositive(limit)) {
            if (isTRUE(perMillion)) {
                if (limit < 1e6L) {
                    warning("`limit`: Use absolute value, not per million.")
                } else {
                    limit <- limit / 1e6L
                }
            }
            p <- p + acid_geom_abline(yintercept = limit)
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

formals(plotTotalReads.bcbioRNASeq)[["fill"]] <-
    formalsList[["fill.discrete"]]
formals(plotTotalReads.bcbioRNASeq)[["flip"]] <-
    formalsList[["flip"]]


#' @rdname plotTotalReads
#' @export
setMethod(
    f = "plotTotalReads",
    signature = signature("bcbioRNASeq"),
    definition = plotTotalReads.bcbioRNASeq
)
