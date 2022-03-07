#' @name plotExonicMappingRate
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inherit AcidGenerics::plotExonicMappingRate
#' @note Updated 2022-03-07.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @description
#' Ideally, at least 60 percent of total reads should map to exons for RNA-seq.
#'
#' @examples
#' data(bcb)
#'
#' ## bcbioRNASeq ====
#' plotExonicMappingRate(bcb)
NULL



## Updated 2022-03-07.
`plotExonicMappingRate,bcbioRNASeq` <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        limit = 0.6,
        labels = list(
            "title" = "Exonic mapping rate",
            "subtitle" = NULL,
            "sampleAxis" = NULL,
            "metricAxis" = "exonic mapping rate (%)"
        ),
        flip = getOption(x = "acid.flip", default = TRUE)
    ) {
        validObject(object)
        assert(
            isProportion(limit),
            isFlag(flip)
        )
        labels <- matchLabels(labels)
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
        labels[["fill"]] <- paste(interestingGroups, collapse = ":\n")
        names(labels)[names(labels) == "sampleAxis"] <- "x"
        names(labels)[names(labels) == "metricAxis"] <- "y"
        p <- p + do.call(what = labs, args = labels)
        ## Limit.
        if (isPositive(limit)) {
            limit <- limit * 100L
            if (limit < 100L) {
                p <- p + acid_geom_abline(yintercept = limit)
            }
        }
        ## Color palette.
        p <- p + autoDiscreteFillScale()
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

formals(`plotExonicMappingRate,bcbioRNASeq`)[["flip"]] <-
    formalsList[["flip"]]



#' @rdname plotExonicMappingRate
#' @export
setMethod(
    f = "plotExonicMappingRate",
    signature = signature(object = "bcbioRNASeq"),
    definition = `plotExonicMappingRate,bcbioRNASeq`
)
