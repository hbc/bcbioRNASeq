#' @name plotGeneSaturation
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inherit AcidGenerics::plotGeneSaturation
#' @note Updated 2021-07-21.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' ## bcbioRNASeq ====
#' data(bcb)
#' plotGeneSaturation(bcb, label = FALSE)
#' plotGeneSaturation(bcb, label = TRUE)
NULL



## Updated 2021-09-10.
`plotGeneSaturation,bcbioRNASeq` <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        minCounts = 1L,
        perMillion = TRUE,
        trendline = FALSE,
<<<<<<< HEAD
<<<<<<< HEAD
        label = getOption(
            x = "acid.label",
            default = FALSE
        ),
=======
        label,
>>>>>>> 58c64d607b3b (Improve color handling)
=======
        label = getOption(x = "acid.label", default = FALSE),
>>>>>>> 19ac90113549 (Simplify code depending on "formalsList")
        labels = list(
            "title" = "Gene saturation",
            "subtitle" = NULL,
            "x" = "mapped reads",
            "y" = "gene count"
        )
    ) {
        validObject(object)
        assert(
            .isGeneLevel(object),
            isInt(minCounts),
            isInRange(minCounts, lower = 1L, upper = Inf),
            isFlag(perMillion),
            isFlag(trendline),
<<<<<<< HEAD
<<<<<<< HEAD
            isFlag(label)
=======
            isFlag(label),
            isGGScale(color, scale = "discrete", aes = "color", nullOK = TRUE)
>>>>>>> 23d4c4234bea (Draft update to label matching)
=======
            isFlag(label)
>>>>>>> 58c64d607b3b (Improve color handling)
        )
        labels <- matchLabels(labels)
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
        labels[["color"]] <- paste(interestingGroups, collapse = ":\n")
        p <- p + do.call(what = labs, args = labels)
        ## Trendline.
        if (isTRUE(trendline)) {
            p <- p + geom_smooth(method = "lm", se = FALSE)
        }
        ## Color palette.
        p <- p + autoDiscreteColorScale()
        ## Hide sample name legend.
        if (isTRUE(label)) {
            p <- p + acid_geom_label_repel(
                mapping = aes(label = !!sym("sampleName"))
            )
        }
        ## Return.
        p
    }

<<<<<<< HEAD
<<<<<<< HEAD
=======
formals(`plotGeneSaturation,bcbioRNASeq`)[["label"]] <-
    formalsList[["label"]]

>>>>>>> 58c64d607b3b (Improve color handling)
=======
>>>>>>> 19ac90113549 (Simplify code depending on "formalsList")


#' @rdname plotGeneSaturation
#' @export
setMethod(
    f = "plotGeneSaturation",
    signature = signature(object = "bcbioRNASeq"),
    definition = `plotGeneSaturation,bcbioRNASeq`
)
