#' @name plotMappedReads
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inherit bioverbs::plotMappedReads
#' @inheritParams minimalism::params
#' @inheritParams basejump::params
#' @inheritParams params
#'
#' @examples
#' data(bcb)
#' plotMappedReads(bcb)
NULL



#' @rdname plotMappedReads
#' @name plotMappedReads
#' @importFrom bioverbs plotMappedReads
#' @export
NULL



plotMappedReads.bcbioRNASeq <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        limit = 10e6L,
        perMillion = TRUE,
        fill,
        flip,
        title = "mapped reads"
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

        # Convert to per million, if desired.
        yLab <- "reads"
        if (isTRUE(perMillion)) {
            data[["mappedReads"]] <- data[["mappedReads"]] / 1e6L
            yLab <- paste(yLab, "per million")
        }

        p <- ggplot(
            data = data,
            mapping = aes(
                x = !!sym("sampleName"),
                y = !!sym("mappedReads"),
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
            p <- flip_x_discrete(p)
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
