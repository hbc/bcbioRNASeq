#' @name plotRRNAMappingRate
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inherit acidgenerics::plotRRNAMappingRate
#' @note Updated 2019-09-16.
#'
#' @inheritParams acidroxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' data(bcb)
#' plotRRNAMappingRate(bcb)
NULL



#' @rdname plotRRNAMappingRate
#' @name plotRRNAMappingRate
#' @importFrom acidgenerics plotRRNAMappingRate
#' @usage plotRRNAMappingRate(object, ...)
#' @export
NULL



## Updated 2019-09-16.
`plotRRNAMappingRate,bcbioRNASeq` <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        limit = 0.1,
        fill,
        labels = list(
            title = "rRNA mapping rate",
            subtitle = NULL,
            sampleAxis = NULL,
            metricAxis = "rRNA mapping rate (%)"
        ),
        flip
    ) {
        validObject(object)
        assert(
            isNumber(limit),
            isProportion(limit),
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
        metricCol <- "rrnaRate"
        ## Warn and early return if rRNA rate was not calculated.
        if (!isSubset(metricCol, colnames(data))) {
            warning("rRNA mapping rate was not calculated. Skipping plot.")
            return(invisible())
        }
        p <- ggplot(
            data = data,
            mapping = aes(
                x = !!sym("sampleName"),
                y = !!sym(metricCol) * 100L,
                fill = !!sym("interestingGroups")
            )
        ) +
            acid_geom_bar() +
            acid_scale_y_continuous_nopad(limits = c(0L, 100L))
        ## Labels.
        if (is.list(labels)) {
            labels[["fill"]] <- paste(interestingGroups, collapse = ":\n")
            names(labels)[names(labels) == "sampleAxis"] <- "x"
            names(labels)[names(labels) == "metricAxis"] <- "y"
            p <- p + do.call(what = labs, args = labels)
        }
        ## Limit.
        if (isPositive(limit)) {
            ## Convert to percentage
            if (limit > 1L) {
                ## nocov start
                warning("'limit': Use ratio (0-1) instead of percentage.")
                ## nocov end
            } else {
                limit <- limit * 100L
            }
            if (limit < 100L) {
                p <- p + acid_geom_abline(yintercept = limit)
            }
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

formals(`plotRRNAMappingRate,bcbioRNASeq`)[c("fill", "flip")] <-
    formalsList[c("fill.discrete", "flip")]



#' @rdname plotRRNAMappingRate
#' @export
setMethod(
    f = "plotRRNAMappingRate",
    signature = signature("bcbioRNASeq"),
    definition = `plotRRNAMappingRate,bcbioRNASeq`
)
