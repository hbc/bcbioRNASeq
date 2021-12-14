#' @name plotTotalReads
#' @author Michael Steinbaugh, Rory Kirchner, Victor Barrera
#' @inherit AcidGenerics::plotTotalReads
#' @note Updated 2021-09-10.
#'
#' @inheritParams AcidRoxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' ## bcbioRNASeq ====
#' data(bcb)
#' plotTotalReads(bcb)
NULL



## Updated 2021-09-10.
`plotTotalReads,bcbioRNASeq` <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        limit = 20e6L,
        perMillion = TRUE,
        labels = list(
            "title" = "Total reads",
            "subtitle" = NULL,
            "sampleAxis" = NULL,
            "metricAxis" = "reads"
        ),
<<<<<<< HEAD
        flip = getOption(
            x = "acid.flip",
            default = TRUE
        )
=======
        flip = getOption(x = "acid.flip", default = TRUE)
>>>>>>> 19ac90113549 (Simplify code depending on "formalsList")
    ) {
        validObject(object)
        assert(
            isInt(limit),
            isNonNegative(limit),
            isFlag(perMillion),
            isFlag(flip)
        )
        labels <- matchLabels(labels)
        interestingGroups(object) <-
            matchInterestingGroups(object, interestingGroups)
        interestingGroups <- interestingGroups(object)
        data <- metrics(object)
        if (isTRUE(perMillion)) {
            data[["totalReads"]] <- data[["totalReads"]] / 1e6L
            if (isString(labels[["metricAxis"]])) {
                labels[["metricAxis"]] <-
                    paste(labels[["metricAxis"]], "(per million)")
            }
        }
        p <- ggplot(
                data = data,
                mapping = aes(
                    x = !!sym("sampleName"),
                    y = !!sym("totalReads"),
                    fill = !!sym("interestingGroups")
                )
            ) +
            acid_geom_bar() +
            acid_scale_y_continuous_nopad()
        ## Labels.
        labels[["fill"]] <- paste(interestingGroups, collapse = ":\n")
        names(labels)[names(labels) == "sampleAxis"] <- "x"
        names(labels)[names(labels) == "metricAxis"] <- "y"
        p <- p + do.call(what = labs, args = labels)
        ## Limit.
        if (isPositive(limit)) {
            if (isTRUE(perMillion)) {
                assert(
                    isTRUE(limit >= 1e6L),
                    msg = sprintf(
                        "'%s': %s",
                        limit,
                        "Use absolute value (1e7), not per million (1)."
                    )
                )
                limit <- limit / 1e6L
            }
            p <- p + acid_geom_abline(yintercept = limit)
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
formals(`plotTotalReads,bcbioRNASeq`)[["flip"]] <-
    formalsList[["flip"]]

>>>>>>> d8f572c77c65 (Finish simplifying color handling in QC plots)
=======
>>>>>>> 19ac90113549 (Simplify code depending on "formalsList")


#' @rdname plotTotalReads
#' @export
setMethod(
    f = "plotTotalReads",
    signature = signature(object = "bcbioRNASeq"),
    definition = `plotTotalReads,bcbioRNASeq`
)
