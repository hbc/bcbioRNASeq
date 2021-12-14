#' @name plotMappingRate
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inherit AcidGenerics::plotMappingRate
#' @note Updated 2021-09-10.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' ## bcbioRNASeq ====
#' data(bcb)
#' plotMappingRate(bcb)
NULL



## Updated 2021-07-21.
`plotMappingRate,bcbioRNASeq` <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        limit = 0.7,
        labels = list(
            "title" = "Mapping rate",
            "subtitle" = NULL,
            "sampleAxis" = NULL,
            "metricAxis" = "mapping rate (%)"
        ),
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> 17bec6711268 (Draft update)
        flip = getOption(
            x = "acid.flip",
            default = TRUE
        )
<<<<<<< HEAD
=======
        flip = getOption(x = "acid.flip", default = TRUE)
>>>>>>> 19ac90113549 (Simplify code depending on "formalsList")
=======
>>>>>>> 17bec6711268 (Draft update)
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
        p <- ggplot(
            data = metrics(object),
            mapping = aes(
                x = !!sym("sampleName"),
                y = !!sym("mappedReads") / !!sym("totalReads") * 100L,
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

<<<<<<< HEAD
<<<<<<< HEAD
=======
formals(`plotMappingRate,bcbioRNASeq`)[["flip"]] <-
    formalsList[["flip"]]

>>>>>>> d8f572c77c65 (Finish simplifying color handling in QC plots)
=======
>>>>>>> 19ac90113549 (Simplify code depending on "formalsList")


#' @rdname plotMappingRate
#' @export
setMethod(
    f = "plotMappingRate",
    signature = signature(object = "bcbioRNASeq"),
    definition = `plotMappingRate,bcbioRNASeq`
)
