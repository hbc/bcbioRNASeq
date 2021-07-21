#' @name plotExonicMappingRate
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inherit AcidGenerics::plotExonicMappingRate
#' @note Updated 2021-07-21.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @description
#' Ideally, at least 60 percent of total reads should map to exons for RNA-seq.
#'
#' @examples
#' data(bcb)
#' plotExonicMappingRate(bcb)
NULL



## Updated 2021-07-21.
`plotExonicMappingRate,bcbioRNASeq` <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        limit = 0.6,
        fill,
        labels = list(
            title = "Exonic mapping rate",
            subtitle = NULL,
            sampleAxis = NULL,
            metricAxis = "exonic mapping rate (%)"
        ),
        flip
    ) {
        validObject(object)
        assert(
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
        p <- ggplot(
            data = data,
            mapping = aes(
                x = !!sym("sampleName"),
                y = !!sym("exonicRate") * 100L,
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
            limit <- limit * 100L
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
            p <- p + guides(fill = "none")
        }
        ## Return.
        p
    }

formals(`plotExonicMappingRate,bcbioRNASeq`)[c("fill", "flip")] <-
    formalsList[c("fill.discrete", "flip")]



#' @rdname plotExonicMappingRate
#' @export
setMethod(
    f = "plotExonicMappingRate",
    signature = signature("bcbioRNASeq"),
    definition = `plotExonicMappingRate,bcbioRNASeq`
)
