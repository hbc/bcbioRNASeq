#' @name plotGeneSaturation
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#'
#' @inherit bioverbs::plotGeneSaturation
#' @inheritParams minimalism::params
#' @inheritParams basejump::params
#' @inheritParams params
#'
#' @param trendline `logical(1)`.
#'   Include a trendline for each group.
#'
#' @examples
#' data(bcb)
#' plotGeneSaturation(bcb, label = FALSE)
#' plotGeneSaturation(bcb, label = TRUE)
NULL



#' @rdname plotGeneSaturation
#' @name plotGeneSaturation
#' @importFrom bioverbs plotGeneSaturation
#' @export
NULL



plotGeneSaturation.bcbioRNASeq <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        minCounts = 1L,
        perMillion = TRUE,
        trendline = FALSE,
        label,
        color,
        title = "Gene saturation"
    ) {
        validObject(object)
        interestingGroups(object) <-
            matchInterestingGroups(object, interestingGroups)
        interestingGroups <- interestingGroups(object)
        assert(
            isInt(minCounts),
            isInRange(minCounts, lower = 1L, upper = Inf),
            isFlag(perMillion),
            isFlag(trendline),
            isFlag(label),
            isGGScale(color, scale = "discrete", aes = "colour", nullOK = TRUE),
            isString(title, nullOK = TRUE)
        )

        counts <- counts(object, normalized = FALSE)
        data <- metrics(object)
        assert(identical(colnames(counts), data[["sampleID"]]))
        data[["geneCount"]] <- colSums(counts >= minCounts)

        # Convert to per million, if desired.
        xLab <- "mapped reads"
        if (isTRUE(perMillion)) {
            data <- mutate(data, mappedReads = !!sym("mappedReads") / 1e6L)
            xLab <- paste(xLab, "per million")
        }

        p <- ggplot(
            data = data,
            mapping = aes(
                x = !!sym("mappedReads"),
                y = !!sym("geneCount"),
                colour = !!sym("interestingGroups")
            )
        ) +
            geom_point(size = 3L) +
            scale_y_continuous(breaks = pretty_breaks()) +
            expand_limits(x = 0L, y = 0L) +
            labs(
                title = title,
                x = xLab,
                y = "gene count",
                colour = paste(interestingGroups, collapse = ":\n")
            )

        if (isTRUE(trendline)) {
            p <- p + geom_smooth(method = "lm", se = FALSE)
        }

        if (is(color, "ScaleDiscrete")) {
            p <- p + color
        }

        if (isTRUE(label)) {
            p <- p + acid_geom_label_repel(
                mapping = aes(label = !!sym("sampleName"))
            )
        }

        p
    }

formals(plotGeneSaturation.bcbioRNASeq)[["color"]] <-
    formalsList[["color.discrete"]]
formals(plotGeneSaturation.bcbioRNASeq)[["label"]] <-
    formalsList[["label"]]



#' @rdname plotGeneSaturation
#' @export
setMethod(
    f = "plotGeneSaturation",
    signature = signature("bcbioRNASeq"),
    definition = plotGeneSaturation.bcbioRNASeq
)
