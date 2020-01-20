#' @name plotMappedReads
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inherit acidgenerics::plotMappedReads
#' @note Updated 2019-09-16.
#'
#' @inheritParams acidroxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' data(bcb)
#' plotMappedReads(bcb)
NULL



#' @rdname plotMappedReads
#' @name plotMappedReads
#' @importFrom acidgenerics plotMappedReads
#' @usage plotMappedReads(object, ...)
#' @export
NULL



## Updated 2019-09-16.
`plotMappedReads,bcbioRNASeq` <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        limit = 10e6L,
        perMillion = TRUE,
        fill,
        labels = list(
            title = "Mapped reads",
            subtitle = NULL,
            sampleAxis = NULL,
            metricAxis = "reads"
        ),
        flip
    ) {
        validObject(object)
        assert(
            isInt(limit),
            isNonNegative(limit),
            isFlag(perMillion),
            isGGScale(fill, scale = "discrete", aes = "fill", nullOK = TRUE),
            isFlag(flip)
        )
        labels <- matchLabels(
            labels = labels,
            choices = eval(formals()[["labels"]])
        )
        interestingGroups(object) <-
            matchInterestingGroups(object, interestingGroups)
        interestingGroups <- interestingGroups(object)
        data <- metrics(object)
        ## Convert to per million, if desired.
        if (isTRUE(perMillion)) {
            data[["mappedReads"]] <- data[["mappedReads"]] / 1e6L
            if (isString(labels[["metricAxis"]])) {
                labels[["metricAxis"]] <-
                    paste(labels[["metricAxis"]], "(per million)")
            }
        }
        p <- ggplot(
            data = data,
            mapping = aes(
                x = !!sym("sampleName"),
                y = !!sym("mappedReads"),
                fill = !!sym("interestingGroups")
            )
        ) +
            acid_geom_bar() +
            acid_scale_y_continuous_nopad()
        ## Labels.
        if (is.list(labels)) {
            labels[["fill"]] <- paste(interestingGroups, collapse = ":\n")
            names(labels)[names(labels) == "sampleAxis"] <- "x"
            names(labels)[names(labels) == "metricAxis"] <- "y"
            p <- p + do.call(what = labs, args = labels)
        }
        ## Limit.
        if (isPositive(limit)) {
            if (isTRUE(perMillion)) {
                if (limit < 1e6L) {
                    warning("'limit': Use absolute value, not per million.")
                } else {
                    limit <- limit / 1e6L
                }
            }
            p <- p + acid_geom_abline(yintercept = limit)
        }
        ## Fill.
        if (is(fill, "ScaleDiscrete")) {
            p <- p + fill
        }
        ## Flip.
        if (isTRUE(flip)) {
            p <- acid_coord_flip(p)
        }
        ## Hide sample name legend.
        if (identical(interestingGroups, "sampleName")) {
            p <- p + guides(fill = FALSE)
        }
        ## Return.
        p
    }

formals(`plotMappedReads,bcbioRNASeq`)[c("fill", "flip")] <-
    formalsList[c("fill.discrete", "flip")]



#' @rdname plotMappedReads
#' @export
setMethod(
    f = "plotMappedReads",
    signature = signature("bcbioRNASeq"),
    definition = `plotMappedReads,bcbioRNASeq`
)
