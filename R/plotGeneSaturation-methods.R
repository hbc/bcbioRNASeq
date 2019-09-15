#' @name plotGeneSaturation
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inherit bioverbs::plotGeneSaturation
#' @note Updated 2019-09-15.
#'
#' @inheritParams acidroxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' data(bcb)
#' plotGeneSaturation(bcb, label = FALSE)
#' plotGeneSaturation(bcb, label = TRUE)
NULL



#' @rdname plotGeneSaturation
#' @name plotGeneSaturation
#' @importFrom bioverbs plotGeneSaturation
#' @usage plotGeneSaturation(object, ...)
#' @export
NULL



## Updated 2019-09-15.
`plotGeneSaturation,bcbioRNASeq` <-  # nolint
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
        assert(
            .isGeneLevel(object),
            isInt(minCounts),
            isInRange(minCounts, lower = 1L, upper = Inf),
            isFlag(perMillion),
            isFlag(trendline),
            isFlag(label),
            isGGScale(color, scale = "discrete", aes = "color", nullOK = TRUE),
            isString(title, nullOK = TRUE)
        )
        interestingGroups(object) <-
            matchInterestingGroups(object, interestingGroups)
        interestingGroups <- interestingGroups(object)
        counts <- counts(object, normalized = FALSE)
        data <- metrics(object)
        assert(identical(colnames(counts), data[["sampleID"]]))
        data[["geneCount"]] <- colSums(counts >= minCounts)
        ## Convert to per million, if desired.
        xLab <- "mapped reads"
        if (isTRUE(perMillion)) {
            data[["mappedReads"]] <- data[["mappedReads"]] / 1e6L
            xLab <- paste(xLab, "per million")
        }
        p <- ggplot(
            data = data,
            mapping = aes(
                x = !!sym("mappedReads"),
                y = !!sym("geneCount"),
                color = !!sym("interestingGroups")
            )
        ) +
            geom_point(size = 3L) +
            scale_y_continuous(breaks = pretty_breaks()) +
            expand_limits(x = 0L, y = 0L) +
            labs(
                title = title,
                x = xLab,
                y = "gene count",
                color = paste(interestingGroups, collapse = ":\n")
            )
        ## Trendline.
        if (isTRUE(trendline)) {
            p <- p + geom_smooth(method = "lm", se = FALSE)
        }
        ## Color.
        if (is(color, "ScaleDiscrete")) {
            p <- p + color
        }
        ## Hide sample name legend.
        if (isTRUE(label)) {
            p <- p + acid_geom_label_repel(
                mapping = aes(label = !!sym("sampleName"))
            )
        }
        ## Return.
        p
    }

formals(`plotGeneSaturation,bcbioRNASeq`)[c("color", "label")] <-
    formalsList[c("color.discrete", "label")]



#' @rdname plotGeneSaturation
#' @export
setMethod(
    f = "plotGeneSaturation",
    signature = signature("bcbioRNASeq"),
    definition = `plotGeneSaturation,bcbioRNASeq`
)
