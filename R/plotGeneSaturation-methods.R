#' @name plotGeneSaturation
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inherit AcidGenerics::plotGeneSaturation
#' @note Updated 2022-05-09.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' data(bcb)
#'
#' ## bcbioRNASeq ====
#' plotGeneSaturation(bcb, label = FALSE)
#' plotGeneSaturation(bcb, label = TRUE)
NULL



## Updated 2022-05-09.
`plotGeneSaturation,bcbioRNASeq` <- # nolint
    function(object,
             interestingGroups = NULL,
             minCounts = 1L,
             perMillion = TRUE,
             trendline = FALSE,
             label = getOption(x = "acid.label", default = FALSE),
             labels = list(
                 "title" = "Gene saturation",
                 "subtitle" = NULL,
                 "x" = "mapped reads",
                 "y" = "gene count"
             )) {
        validObject(object)
        assert(
            .isGeneLevel(object),
            isInt(minCounts),
            isInRange(minCounts, lower = 1L, upper = Inf),
            isFlag(perMillion),
            isFlag(trendline),
            isFlag(label)
        )
        labels <- matchLabels(labels)
        interestingGroups(object) <-
            matchInterestingGroups(object, interestingGroups)
        interestingGroups <- interestingGroups(object)
        counts <- counts(object, normalized = FALSE)
        data <- metrics(object)
        assert(identical(colnames(counts), rownames(data)))
        data[["geneCount"]] <- colSums(counts >= minCounts)
        ## Convert to per million, if desired.
        if (isTRUE(perMillion)) {
            data[["mappedReads"]] <- data[["mappedReads"]] / 1e6L
            if (isString(labels[["x"]])) {
                labels[["x"]] <- paste(labels[["x"]], "(per million)")
            }
        }
        p <- ggplot(
            data = as.data.frame(data),
            mapping = aes(
                x = .data[["mappedReads"]],
                y = .data[["geneCount"]],
                color = .data[["interestingGroups"]]
            )
        ) +
            geom_point(size = 3L) +
            scale_y_continuous(breaks = pretty_breaks()) +
            expand_limits(x = 0L, y = 0L)
        ## Labels.
        labels[["color"]] <- paste(interestingGroups, collapse = ":\n")
        p <- p + do.call(what = labs, args = labels)
        ## Trendline.
        if (isTRUE(trendline)) {
            p <- p + geom_smooth(method = "lm", se = FALSE)
        }
        ## Color palette.
        p <- p + acid_scale_color_discrete()
        ## Hide sample name legend.
        if (isTRUE(label)) {
            p <- p + acid_geom_label_repel(
                mapping = aes(label = .data[["sampleName"]])
            )
        }
        ## Return.
        p
    }



#' @rdname plotGeneSaturation
#' @export
setMethod(
    f = "plotGeneSaturation",
    signature = signature(object = "bcbioRNASeq"),
    definition = `plotGeneSaturation,bcbioRNASeq`
)
