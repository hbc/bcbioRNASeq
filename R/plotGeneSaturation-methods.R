#' @name plotGeneSaturation
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inherit AcidGenerics::plotGeneSaturation
#' @note Updated 2019-09-16.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' data(bcb)
#' plotGeneSaturation(bcb, label = FALSE)
#' plotGeneSaturation(bcb, label = TRUE)
NULL



## Updated 2019-09-16.
`plotGeneSaturation,bcbioRNASeq` <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        minCounts = 1L,
        perMillion = TRUE,
        trendline = FALSE,
        label,
        color,
        labels = list(
            title = "Gene saturation",
            subtitle = NULL,
            x = "mapped reads",
            y = "gene count"
        )
    ) {
        validObject(object)
        assert(
            .isGeneLevel(object),
            isInt(minCounts),
            isInRange(minCounts, lower = 1L, upper = Inf),
            isFlag(perMillion),
            isFlag(trendline),
            isFlag(label),
            isGGScale(color, scale = "discrete", aes = "color", nullOK = TRUE)
        )
        labels <- matchLabels(
            labels = labels,
            choices = eval(formals()[["labels"]])
        )
        interestingGroups(object) <-
            matchInterestingGroups(object, interestingGroups)
        interestingGroups <- interestingGroups(object)
        counts <- counts(object, normalized = FALSE)
        data <- metrics(object)
        assert(identical(colnames(counts), data[["sampleId"]]))
        data[["geneCount"]] <- colSums(counts >= minCounts)
        ## Convert to per million, if desired.
        if (isTRUE(perMillion)) {
            data[["mappedReads"]] <- data[["mappedReads"]] / 1e6L
            if (isString(labels[["x"]])) {
                labels[["x"]] <- paste(labels[["x"]], "(per million)")
            }
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
            expand_limits(x = 0L, y = 0L)
        ## Labels.
        if (is.list(labels)) {
            labels[["color"]] <- paste(interestingGroups, collapse = ":\n")
            p <- p + do.call(what = labs, args = labels)
        }
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
