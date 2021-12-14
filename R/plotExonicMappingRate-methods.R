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
#' ## bcbioRNASeq ====
#' data(bcb)
#' plotExonicMappingRate(bcb)
NULL



## Updated 2021-09-10.
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
<<<<<<< HEAD
<<<<<<< HEAD
        ## Flip.
=======
        ## Flip, if desired.
>>>>>>> 58c64d607b3b (Improve color handling)
=======
        ## Flip.
>>>>>>> 88e8e9be9797 (Improve default color handling)
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
formals(`plotExonicMappingRate,bcbioRNASeq`)[["flip"]] <-
    formalsList[["flip"]]

>>>>>>> 58c64d607b3b (Improve color handling)
=======
>>>>>>> 19ac90113549 (Simplify code depending on "formalsList")


#' @rdname plotExonicMappingRate
#' @export
setMethod(
    f = "plotExonicMappingRate",
    signature = signature(object = "bcbioRNASeq"),
    definition = `plotExonicMappingRate,bcbioRNASeq`
)
